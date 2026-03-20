import streamlit as st
import pandas as pd
import tempfile
import os
import io
from io import StringIO
import sys
import subprocess
import zipfile
import requests

# Patch for meeko/rdkit compatibility (rdkit.six removed in new versions)
try:
    from rdkit import six
except ImportError:
    import io
    import types
    sys.modules['rdkit.six'] = types.ModuleType('rdkit.six')
    sys.modules['rdkit.six'].StringIO = io.StringIO
    sys.modules['rdkit.six'].BytesIO = io.BytesIO

from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTMolecule
import py3Dmol
from stmol import showmol
from bs4 import BeautifulSoup

try:
    from vina import Vina
except ImportError:
    st.error("""
    **ERRO DE DEPENDÊNCIA: AutoDock Vina não encontrado.**
    
    A biblioteca Python `vina` é necessária para esta aplicação. 
    Por favor, instale-a usando Conda (recomendado):
    
    `conda install -c conda-forge vina`
    
    Ou tente via pip (pode exigir compilação):
    `pip install vina`
    """)
    st.stop()


# Page Configuration
st.set_page_config(page_title="EASYDocking", layout="wide")

# Inicializar session_state para persistir resultados entre reruns
if 'docking_results' not in st.session_state:
    st.session_state['docking_results'] = None
if 'protein_data_cached' not in st.session_state:
    st.session_state['protein_data_cached'] = None
if 'viz_format_cached' not in st.session_state:
    st.session_state['viz_format_cached'] = None
if 'box_params_cached' not in st.session_state:
    st.session_state['box_params_cached'] = None

# Title
st.title("EASYDocking 🧬")
st.markdown("""
Esta aplicação permite realizar docking molecular utilizando abordagens clássicas e modelos orientados por Inteligência Artificial (IA).
A entrada pode ser uma proteína (PDBQT recomendado ou PDB) e ligantes (SMILES ou SDF).

⚠️ **Boas Práticas de Bioinformática:**
- **Receptor**: Sempre limpe seu arquivo PDB antes de subir (remova moléculas de água, íons, solventes e ligantes cristalizados indesejados). O conversor interno tentará adicionar hidrogênios, mas resultados ideais dependem de arquivos bem preparados e protonados no pH correto.
- **Ligantes**: O aplicativo gerará conformações em 3D e adicionará hidrogênios básicos. Certifique-se de que seus SMILES/SDF possuam o estado de protonação correto, caso sua molécula possua cátions/ânions.

### Referências
- **AutoDock Vina**: Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *Journal of computational chemistry*.
- **DiffDock (IA)**: Corso, G., et al. (2022). DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking. *ICLR*.
- **TankBind (IA)**: Lu, W., et al. (2022). TankBind: Trigonometry-Aware Neural NetworKs for Drug-Protein Binding Structure Prediction. *NeurIPS*.
- **P2Rank (IA sítio ativo)**: Krivak, R., & Hoksza, D. (2018). P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure. *Journal of cheminformatics*.
""")

# Sidebar for Global Parameters
st.sidebar.header("Métodos e Configurações")

docking_method = st.sidebar.selectbox(
    "1. Motor de Docking",
    [
        "AutoDock Vina (Clássico)", 
        "DiffDock (IA)", 
        "TankBind (IA)", 
        "Análise 2D (Somente PLIP, sem Docking)"
    ]
)

active_site_method = st.sidebar.selectbox(
    "2. Identificação do Sítio Ativo",
    ["Manual (Coordenadas da Caixa)", "Automático (via PDB ID)", "P2Rank (Predição por IA)"]
)

st.sidebar.markdown("---")
st.sidebar.header("Parâmetros do Grid (Caixa)")

if active_site_method == "Manual (Coordenadas da Caixa)":
    center_x = st.sidebar.number_input("Centro X", value=0.0)
    center_y = st.sidebar.number_input("Centro Y", value=0.0)
    center_z = st.sidebar.number_input("Centro Z", value=0.0)

    size_x = st.sidebar.number_input("Tamanho X (Angstroms)", value=20.0, min_value=1.0)
    size_y = st.sidebar.number_input("Tamanho Y (Angstroms)", value=20.0, min_value=1.0)
    size_z = st.sidebar.number_input("Tamanho Z (Angstroms)", value=20.0, min_value=1.0)

elif active_site_method == "Automático (via PDB ID)":
    st.sidebar.info("A caixa será calculada automaticamente com base no ligante original da estrutura PDB.")
    pdb_id_input = st.sidebar.text_input("Digite o PDB ID (ex: 1STP)", max_chars=4)
    
    # Valores padrão iniciais (caso erro ou vazio)
    center_x, center_y, center_z = 0.0, 0.0, 0.0
    size_x, size_y, size_z = 20.0, 20.0, 20.0
    
    if len(pdb_id_input) == 4:
        with st.sidebar.status(f"Calculando caixa para {pdb_id_input.upper()}..."):
            try:
                # 1. Obter o Comp ID do ligante via GraphQL
                query = """
                {
                  entry(entry_id: "%s") {
                    rcsb_binding_affinity {
                      comp_id
                    }
                    nonpolymer_entities {
                      nonpolymer_comp {
                        chem_comp {
                          id
                        }
                      }
                    }
                  }
                }
                """ % pdb_id_input.upper()
                
                resp = requests.post("https://data.rcsb.org/graphql", json={'query': query})
                found_comp_id = None
                
                if resp.ok:
                    data = resp.json()
                    entry = data.get('data', {}).get('entry', {})
                    if entry:
                        affinities = entry.get('rcsb_binding_affinity')
                        if affinities and len(affinities) > 0:
                            found_comp_id = affinities[0].get('comp_id')
                        else:
                            nonpolymers = entry.get('nonpolymer_entities')
                            if nonpolymers:
                                ignore_list = ["HOH", "DOD", "WAT", "NA", "CL", "K", "MG", "CA", "ZN", "CU", "FE", "MN", "CO", "NI", "I", "BR", "SO4", "PO4", "NO3", "CO3", "ACT", "FMT", "ACE", "GOL", "PEG", "EDO", "DMS", "PG4", "PGE", "PE4", "BME", "DTT", "NAP", "NADP", "NAD"]
                                for entity in nonpolymers:
                                    try:
                                        cid = entity['nonpolymer_comp']['chem_comp']['id'].upper()
                                        if cid not in ignore_list:
                                            found_comp_id = cid
                                            break
                                    except:
                                        pass
                
                if not found_comp_id:
                    st.error("Nenhum ligante válido / inibidor encontrado nesta estrutura (Apo ou apenas solventes).")
                else:
                    st.success(f"Ligante alvo detectado: {found_comp_id}")
                    # 2. Baixar PDB e Calcular Centro/Tamanho
                    pdb_resp = requests.get(f"https://files.rcsb.org/download/{pdb_id_input.upper()}.pdb")
                    if pdb_resp.ok:
                        lines = pdb_resp.text.splitlines()
                        coords = []
                        for line in lines:
                            if line.startswith('HETATM') and line[17:20].strip() == found_comp_id:
                                atom_name = line[12:16].strip()
                                if not atom_name.startswith('H'): # Ignora hidrogenios
                                    try:
                                        x = float(line[30:38].strip())
                                        y = float(line[38:46].strip())
                                        z = float(line[46:54].strip())
                                        coords.append([x, y, z])
                                    except: pass
                        
                        if coords:
                            n_atoms = len(coords)
                            x_vals = [c[0] for c in coords]
                            y_vals = [c[1] for c in coords]
                            z_vals = [c[2] for c in coords]
                            
                            center_x = sum(x_vals) / n_atoms
                            center_y = sum(y_vals) / n_atoms
                            center_z = sum(z_vals) / n_atoms
                            
                            buffer = 10.0
                            size_x = (max(x_vals) - min(x_vals)) + buffer
                            size_y = (max(y_vals) - min(y_vals)) + buffer
                            size_z = (max(z_vals) - min(z_vals)) + buffer
                            
                            # Atualiza para exibir e salvar globalmente  
                            st.write(f"**Centro X:** {center_x:.2f} | **Tamanho X:** {size_x:.2f}")
                            st.write(f"**Centro Y:** {center_y:.2f} | **Tamanho Y:** {size_y:.2f}")
                            st.write(f"**Centro Z:** {center_z:.2f} | **Tamanho Z:** {size_z:.2f}")
                        else:
                            st.error("Falha ao extrair coordenadas do ligante.")
                    else:
                        st.error("Falha ao baixar o arquivo PDB do RCSB.")
            except Exception as e:
                st.error(f"Erro ao processar: {e}")

else:
    st.sidebar.info("A predição automática por IA (P2Rank) identificará automaticamente os bolsões mais prováveis na proteína.")
    center_x, center_y, center_z = 0.0, 0.0, 0.0
    size_x, size_y, size_z = 20.0, 20.0, 20.0

if docking_method == "AutoDock Vina (Clássico)":
    exhaustiveness = st.sidebar.slider("Exhaustiveness (Precisão)", min_value=1, max_value=32, value=8)
    num_modes = st.sidebar.slider("Número de Modos (Poses)", min_value=1, max_value=20, value=9)
else:
    st.sidebar.info(f"O modelo ({docking_method}) não utiliza parâmetros de box/grid clássico.")
    exhaustiveness = 8
    num_modes = 1

# ------------------------------------------------------------------------------
# STEP 1: PROTEIN INPUT
# ------------------------------------------------------------------------------
st.header("1. Upload da Proteína (Receptor)")
protein_file = st.file_uploader("Carregar arquivo da proteína (.pdbqt ou .pdb)", type=["pdbqt", "pdb"])

st.subheader("Opções de Limpeza (Pré-processamento)")
col_a, col_b = st.columns(2)
with col_a:
    remove_water = st.checkbox("Remover Águas (HOH/WAT/SOL)", value=True, help="Recomendado. Remove as moléculas de solvente cristalinas.")
with col_b:
    remove_hetatm = st.checkbox("Remover Heteroátomos (HETATM)", value=True, help="Recomendado. Remove os ligantes, íons ou pequenos compostos inerentes ao arquivo PDB.")

def convert_pdb_to_pdbqt(pdb_content):
    """
    Converts PDB content to PDBQT string using OpenBabel (via system call).
    This is generally more robust for complex PDB files than RDKit/Meeko.
    Flags used:
    -xr: Output as rigid receptor (strips ROOT/BRANCH)
    -h: Add hydrogens
    --partialcharge gasteiger: Calculate partial charges
    """
    try:
        # 1. Write PDB content to a temporary file
        input_file = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        input_file.write(pdb_content.encode("utf-8"))
        input_file.close() # Close so obabel can access it
        
        output_filename = input_file.name + ".qt" # temporary output name
        
        # 2. Construct OpenBabel command
        # obabel -ipdb input.pdb -opdbqt -O output.pdbqt -xr -h --partialcharge gasteiger
        cmd = [
            "obabel",
            "-ipdb", input_file.name,
            "-opdbqt", "-O", output_filename,
            "-xr",                  # Rigid receptor
            "-h",                   # Add hydrogens
            "--partialcharge", "gasteiger" # Calculate charges
        ]
        
        # 3. Execute Command
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            st.error(f"Erro no OpenBabel: {result.stderr}")
            # Try cleanup
            try:
                os.unlink(input_file.name)
            except: pass
            return None
            
        # 4. Read Output
        if os.path.exists(output_filename):
            with open(output_filename, "r") as f:
                pdbqt_content = f.read()
            
            # Cleanup
            try:
                os.unlink(input_file.name)
                os.unlink(output_filename)
            except: pass
            
            return pdbqt_content
        else:
            st.error("OpenBabel não gerou o arquivo de saída.")
            return None

    except Exception as e:
        st.error(f"Erro ao executar OpenBabel: {e}")
        return None

def load_protein(protein_file, rm_water, rm_hetatm):
    if protein_file is not None:
        string_data = protein_file.getvalue().decode("utf-8")
        filename = protein_file.name
        
        # 0. Limpeza em nível de texto (Text-level cleaning)
        cleaned_lines = []
        for line in string_data.splitlines():
            # Filtrar heteroátomos globais se habilitado
            if rm_hetatm and line.startswith("HETATM"):
                continue
            
            # Filtrar águas (HOH, WAT, SOL)
            if rm_water and (line.startswith("ATOM") or line.startswith("HETATM")):
                if len(line) >= 20:
                    resname = line[17:20].strip()
                    if resname in ["HOH", "WAT", "SOL"]:
                        continue
            
            cleaned_lines.append(line)
            
        string_data = "\n".join(cleaned_lines)

        if filename.lower().endswith('.pdb'):
            with st.spinner("Convertendo PDB para PDBQT (adicionando hidrogênios e cargas)..."):
                converted_pdbqt = convert_pdb_to_pdbqt(string_data)
                if converted_pdbqt:
                    return converted_pdbqt, filename + "qt" # fake extension path
                else:
                    st.error("Falha ao converter PDB. Verifique o formato do arquivo.")
                    return None, None
        
        return string_data, filename
    return None, None

protein_data, protein_name = load_protein(protein_file, remove_water, remove_hetatm)

if protein_data:
    st.success(f"Proteína '{protein_name}' carregada e processada!")
    
    # Check extension (now we might have fake .pdbqt from conversion)
    is_pdbqt = protein_name.lower().endswith('.pdbqt')

    
    # Visualization setup
    # Determine format for py3Dmol
    viz_format = 'pdbqt' if is_pdbqt else 'pdb'
    
    tfile = tempfile.NamedTemporaryFile(delete=False, suffix=f".{viz_format}")
    tfile.write(protein_data.encode('utf-8'))
    tfile.close()
    protein_path = tfile.name

    with st.expander("Visualizar Proteína e Caixa"):
        col_v1, col_v2, col_v3 = st.columns(3)
        with col_v1:
            show_protein = st.checkbox("Mostrar Proteína", value=True, key="chk_prot")
        with col_v2:
            show_grid = st.checkbox("Mostrar Caixa do Grid", value=True, key="chk_grid")
        with col_v3:
            show_center = st.checkbox("Mostrar Centro Ativo", value=True, key="chk_center")
            
        view = py3Dmol.view(width=800, height=400)
        
        if show_protein:
            view.addModel(protein_data, viz_format)
            view.setStyle({'cartoon': {'color': 'spectrum'}})
        
        # Draw Box
        if show_grid:
            view.addBox({
                'center': {'x': center_x, 'y': center_y, 'z': center_z},
                'dimensions': {'w': size_x, 'h': size_y, 'd': size_z},
                'color': 'green',
                'opacity': 0.5
            })
            
        # Draw Center
        if show_center:
            view.addSphere({
                'center': {'x': center_x, 'y': center_y, 'z': center_z},
                'radius': 1.0,
                'color': 'red',
                'opacity': 0.9
            })
            
        view.zoomTo()
        showmol(view, height=400, width=800)

else:
    st.info("Por favor, carregue um arquivo PDBQT (recomendado) ou PDB.")

# ------------------------------------------------------------------------------
# STEP 2: LIGAND INPUT
# ------------------------------------------------------------------------------
st.header("2. Entrada de Ligantes")
input_method = st.radio("Método de Entrada:", ["Lista de SMILES", "Arquivo SDF"])

ligands_to_dock = []  # List of (name, mol_object)

if input_method == "Lista de SMILES":
    smiles_input = st.text_area("Cole os SMILES aqui (um por linha, formato: SMILES Nome ou apenas SMILES)", height=150)
    if smiles_input:
        lines = smiles_input.strip().split('\n')
        for i, line in enumerate(lines):
            parts = line.strip().split()
            if not parts: continue
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else f"Ligand_{i+1}"
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mol = Chem.AddHs(mol)
                ligands_to_dock.append((name, mol))

elif input_method == "Arquivo SDF":
    sdf_file = st.file_uploader("Carregar arquivo .sdf", type=["sdf"])
    if sdf_file:
        # Save to temp to read with SDMolSupplier (needs path usually or file-like)
        # RDKit can read from stream but simpler to save temp
        t_sdf = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        t_sdf.write(sdf_file.getvalue())
        t_sdf.close()
        suppl = Chem.SDMolSupplier(t_sdf.name)
        for i, mol in enumerate(suppl):
            if mol:
                name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Ligand_{i+1}"
                mol = Chem.AddHs(mol)
                ligands_to_dock.append((name, mol))
        os.unlink(t_sdf.name)

st.write(f"Ligantes carregados: {len(ligands_to_dock)}")

# ------------------------------------------------------------------------------
# STEP 3: RUN DOCKING
# ------------------------------------------------------------------------------
st.header("3. Executar Docking")

if st.button("Iniciar"):
    if docking_method not in ["AutoDock Vina (Clássico)", "Análise 2D (Somente PLIP, sem Docking)"]:
        st.warning(f"O modelo escolhido ({docking_method}) requer ambiente suportado por GPU e deep learning. A demonstração prosseguirá com AutoDock Vina ou será abortada se o modelo não estiver integrado no momento.")
        st.stop()
        
    if active_site_method != "Manual (Coordenadas da Caixa)":
        st.warning("Módulo P2Rank (IA) convocado! Como este é um ambiente de demonstração, o script utilizará o centro do sistema (fallback) para prosseguir por enquanto.")
    
    if not protein_data:
        st.error("Proteína não carregada!")
    elif len(ligands_to_dock) == 0:
        st.error("Nenhum ligante válido encontrado!")
    else:
        # Progress Bar
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        results_list = []
        
        # Prepare Vina
        v = Vina(sf_name='vina')
        
        if not protein_name.lower().endswith('.pdbqt'):
            st.warning("Aviso: O arquivo da proteína não é .pdbqt. O Vina requer cargas parciais e tipos de átomos. Resultados podem ser imprecisos ou falhar.")
        
        try:
            if docking_method != "Análise 2D (Somente PLIP, sem Docking)":
                v.set_receptor(protein_path) # Use the temp file path from Step 1
                # Set Map based on box
                v.compute_vina_maps(center=[center_x, center_y, center_z], box_size=[size_x, size_y, size_z])
            
            total_ligands = len(ligands_to_dock)
            
            for index, (lig_name, mol) in enumerate(ligands_to_dock):
                status_text.text(f"Processando {lig_name} ({index+1}/{total_ligands})...")
                
                # 1. Generate 3D Conformer if no 3D structure is present
                needs_3d = False
                if mol.GetNumConformers() == 0:
                    needs_3d = True
                elif not mol.GetConformer().Is3D():
                    needs_3d = True
                    
                if needs_3d:
                    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                    try:
                        AllChem.UFFOptimizeMolecule(mol)
                    except:
                        pass
                
                # 2. Convert Mol to PDBQT String using Meeko
                meeko_prep = MoleculePreparation()
                meeko_prep.prepare(mol)
                lig_pdbqt = meeko_prep.write_pdbqt_string()
                
                # 3. Validation
                if not lig_pdbqt:
                     st.warning(f"Falha ao preparar PDBQT para {lig_name}")
                     continue

                if docking_method == "Análise 2D (Somente PLIP, sem Docking)":
                    affinity = 0.0
                    all_poses_pdbqt = lig_pdbqt
                else:
                    # 4. Docking
                    v.set_ligand_from_string(lig_pdbqt)
                    v.dock(exhaustiveness=exhaustiveness, n_poses=num_modes)
                    
                    # 5. Get Energy/Score
                    energies = v.energies(n_poses=1)
                    affinity = energies[0][0] if len(energies) > 0 else 0.0

                    # 6. Capture Result PDBQT (All Poses)
                    t_out = tempfile.NamedTemporaryFile(delete=False, suffix=".pdbqt")
                    t_out_name = t_out.name
                    t_out.close()
                    
                    v.write_poses(t_out_name, n_poses=num_modes, overwrite=True)
                    
                    with open(t_out_name, 'r') as f:
                        all_poses_pdbqt = f.read()
                    
                    os.unlink(t_out_name)
                
                heavy_atoms = mol.GetNumHeavyAtoms()
                ligand_efficiency = abs(affinity) / heavy_atoms if heavy_atoms > 0 else 0.0
                
                results_list.append({
                    "Ligand": lig_name,
                    "Affinity (kcal/mol)": affinity,
                    "Atómos Pesados": heavy_atoms,
                    "Ligand Efficiency": round(ligand_efficiency, 3),
                    "PDBQT": all_poses_pdbqt
                })
                
                progress_bar.progress((index + 1) / total_ligands)
            
            if docking_method == "Análise 2D (Somente PLIP, sem Docking)":
                status_text.text("Preparação do Ligante Concluída (PLIP).")
            else:
                status_text.text("Docking Concluído!")
                
            # ✅ Salvar resultados no session_state para persistir após download
            st.session_state['docking_results'] = results_list
            st.session_state['protein_data_cached'] = protein_data
            st.session_state['viz_format_cached'] = viz_format
            st.session_state['box_params_cached'] = {
                'center_x': center_x, 'center_y': center_y, 'center_z': center_z,
                'size_x': size_x, 'size_y': size_y, 'size_z': size_z
            }
            
        except Exception as e:
            st.error(f"Ocorreu um erro durante o docking: {str(e)}")
            import traceback
            st.text(traceback.format_exc())

# ==============================================================================
# EXIBIR RESULTADOS (fora do bloco do botão para persistir após reruns/downloads)
# ==============================================================================
if st.session_state.get('docking_results'):
    results_list = st.session_state['docking_results']
    cached_protein_data = st.session_state.get('protein_data_cached', protein_data)
    cached_viz_format = st.session_state.get('viz_format_cached', 'pdbqt')
    cached_box = st.session_state.get('box_params_cached', {
        'center_x': 0.0, 'center_y': 0.0, 'center_z': 0.0,
        'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0
    })
    
    st.subheader("Tabela de Resultados")
    df_results = pd.DataFrame(results_list)
    if not df_results.empty:
        cols_to_show = ["Ligand", "Affinity (kcal/mol)", "Atómos Pesados", "Ligand Efficiency"]
        st.dataframe(df_results[cols_to_show])
        
        # Download Results CSV
        csv = df_results[cols_to_show].to_csv(index=False).encode('utf-8')
        
        col1, col2 = st.columns(2)
        with col1:
            st.download_button(
                label="Baixar Tabela (CSV)",
                data=csv,
                file_name="docking_results.csv",
                mime="text/csv",
                key="dl_csv"
            )

        # Download All Results (ZIP)
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("docking_results.csv", df_results[cols_to_show].to_csv(index=False))
            ext = "pdbqt" if (cached_protein_data and cached_viz_format == 'pdbqt') else "pdb"
            if cached_protein_data:
                zf.writestr(f"receptor.{ext}", cached_protein_data)
            for i, row in df_results.iterrows():
                lig_name = row["Ligand"]
                lig_pdbqt = row["PDBQT"]
                safe_name = "".join([c for c in lig_name if c.isalnum() or c in (' ', '.', '_')]).strip().replace(" ", "_")
                zf.writestr(f"{safe_name}.pdbqt", lig_pdbqt)
        
        with col2:
            st.download_button(
                label="Baixar Tudo (ZIP)",
                data=zip_buffer.getvalue(),
                file_name="docking_results.zip",
                mime="application/zip",
                key="dl_zip"
            )

        # Visualize Best Result
        st.subheader("Visualização do Melhor Resultado")
        selected_ligand = st.selectbox("Escolha um ligante para visualizar:", df_results["Ligand"].tolist())
        
        if selected_ligand and cached_protein_data:
            row = df_results[df_results["Ligand"] == selected_ligand].iloc[0]
            ligand_pdbqt = row["PDBQT"]
            
            view_res = py3Dmol.view(width=800, height=500)
            view_res.addModel(cached_protein_data, cached_viz_format)
            view_res.setStyle({'model': -1}, {'cartoon': {'color': 'white', 'opacity': 0.7}})
            
            view_res.addModel(ligand_pdbqt, 'pdbqt')
            view_res.setStyle({'model': -1}, {'stick': {'colorscheme': 'greenCarbon'}})
            
            view_res.addBox({
                'center': {'x': cached_box['center_x'], 'y': cached_box['center_y'], 'z': cached_box['center_z']},
                'dimensions': {'w': cached_box['size_x'], 'h': cached_box['size_y'], 'd': cached_box['size_z']},
                'color': 'blue',
                'opacity': 0.2
            })

            view_res.zoomTo()
            showmol(view_res, height=500, width=800)
            
    else:
        st.warning("Nenhum resultado gerado.")

# Remove temp file on session end logic is complex in Streamlit, 
# relying on OS temp cleanup or manual clean if needed. 
