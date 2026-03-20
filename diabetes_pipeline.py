import streamlit as st
import pandas as pd
import tempfile
import os
import io
import sys
import subprocess
import urllib.request
import zipfile

# Patch for meeko/rdkit compatibility
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
import py3Dmol
from stmol import showmol

try:
    from vina import Vina
except ImportError:
    st.error("ERRO: AutoDock Vina não encontrado.")
    st.stop()

# Configuração da Página
st.set_page_config(page_title="Diabetes Pipeline (DPTD)", layout="wide", page_icon="🩸")

st.title("🩸 Diabetes Protein Target Pipeline")
st.markdown("""
Bem-vindo ao **Diabetes Pipeline**, inspirado no *Diabetes Protein Target Database (DPTD)*. 
Aqui você pode selecionar alvos terapêuticos clássicos para a Diabetes Tipo 2, baixar suas estruturas diretamente do servidor PDB, validá-las e iniciar uma varredura de *Virtual Screening* ou Docking com potenciais ligantes.
""")

# Base de Dados Inspirada no DPTD e Atualizada com Dimensões
diabetes_targets = {
    "DPP-4 (Dipeptidyl Peptidase-4)": {"pdb": "2RGU", "res": 1.7, "class": "Enzyme Inhibition", "desc": "Serine protease that rapidly degrades incretins.", "cx": 68.597, "cy": 68.296, "cz": 73.407, "sx": 20.0},
    "Monooxidase (Kynurenine 3-Monooxygenase)": {"pdb": "5X68", "res": 2.1, "class": "Metabolism", "desc": "Regulator linking inflammation to metabolic dysfunction.", "cx": 31.473, "cy": -11.920, "cz": -60.275, "sx": 20.0},
    "IDE (Insulin-Degrading Enzyme)": {"pdb": "4IFH", "res": 2.2, "class": "Enzyme Inhibition", "desc": "Key enzyme for insulin clearance.", "cx": 99.372, "cy": -0.109, "cz": 42.522, "sx": 20.0},
    "IR (Insulin Receptor)": {"pdb": "4IBM", "res": 3.8, "class": "Receptor Modulation", "desc": "Critical for insulin signaling pathway.", "cx": 4.689, "cy": -8.818, "cz": 6.074, "sx": 20.0},
    "GLUT1 (Glucose Transporter 1)": {"pdb": "5EQG", "res": 3.5, "class": "Glucose Transport", "desc": "Facilitates the transport of glucose across plasma membranes.", "cx": 582.427, "cy": -27.075, "cz": 280.935, "sx": 20.0},
    "PPARγ (PPAR Gamma)": {"pdb": "2PRG", "res": 2.3, "class": "Receptor Modulation", "desc": "Nuclear receptor, key regulator of adipocyte differentiation.", "cx": 60.511, "cy": -3.002, "cz": 39.875, "sx": 20.0},
    "Glucokinase (GCK)": {"pdb": "1V4S", "res": 1.9, "class": "Metabolism", "desc": "Sensor of blood glucose in pancreatic beta cells.", "cx": 39.663, "cy": 16.227, "cz": 62.887, "sx": 20.0},
    "GSK-3β (GSK-3 Beta)": {"pdb": "1Q3W", "res": 2.0, "class": "Metabolism", "desc": "Kinase involved in glycogen metabolism and insulin signaling.", "cx": 21.955, "cy": -19.223, "cz": 8.387, "sx": 20.0},
    "FBP1 (Fructose-1,6-Bisphosphatase)": {"pdb": "1FTA", "res": 2.3, "class": "Metabolism", "desc": "Regulatory enzyme in gluconeogenesis.", "cx": 5.366, "cy": 62.848, "cz": 32.856, "sx": 20.0}
}

st.sidebar.header("🎯 Seleção do Alvo (DPTD)")
selected_target_name = st.sidebar.selectbox("Escolha um alvo da Diabetes:", list(diabetes_targets.keys()))
target_info = diabetes_targets[selected_target_name]

st.sidebar.markdown(f"""
**PDB ID:** `{target_info['pdb']}`  
**Resolução:** `{target_info['res']} Å`  
**Classe:** `{target_info['class']}`  
""")
st.sidebar.info(target_info['desc'])

# Grid Setup Customizado no Sidebar
st.sidebar.header("📐 Configurações da Caixa (Grid)")
center_x = st.sidebar.number_input("Center X", value=target_info.get('cx', 0.0), format="%.3f")
center_y = st.sidebar.number_input("Center Y", value=target_info.get('cy', 0.0), format="%.3f")
center_z = st.sidebar.number_input("Center Z", value=target_info.get('cz', 0.0), format="%.3f")
size_x = st.sidebar.number_input("Size X", value=target_info.get('sx', 20.0), min_value=1.0)
size_y = st.sidebar.number_input("Size Y", value=target_info.get('sx', 20.0), min_value=1.0)
size_z = st.sidebar.number_input("Size Z", value=target_info.get('sx', 20.0), min_value=1.0)
exhaustiveness = st.sidebar.slider("Exhaustiveness", min_value=1, max_value=32, value=8)

# Função para baixar do RCSB
@st.cache_data
def fetch_pdb_from_rcsb(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = urllib.request.urlopen(url)
        return response.read().decode('utf-8')
    except Exception as e:
        return None

def convert_pdb_to_pdbqt(pdb_content):
    try:
        input_file = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        input_file.write(pdb_content.encode("utf-8"))
        input_file.close() 
        output_filename = input_file.name + ".qt"
        cmd = ["obabel", "-ipdb", input_file.name, "-opdbqt", "-O", output_filename, "-xr", "-h", "--partialcharge", "gasteiger"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0 and os.path.exists(output_filename):
            with open(output_filename, "r") as f:
                content = f.read()
            os.unlink(input_file.name)
            os.unlink(output_filename)
            return content
    except:
        pass
    return None

# --- Fase 1: Proteína ---
st.header(f"1. Preparação da Proteína: {selected_target_name}")

if "protein_data" not in st.session_state or st.session_state.get("current_pdb") != target_info['pdb']:
    with st.spinner(f"Baixando {target_info['pdb']} do RCSB PDB e convertendo/protonando via OpenBabel..."):
        raw_pdb = fetch_pdb_from_rcsb(target_info['pdb'])
        if raw_pdb:
            pdbqt_data = convert_pdb_to_pdbqt(raw_pdb)
            if pdbqt_data:
                st.session_state["protein_data"] = pdbqt_data
                st.session_state["raw_pdb"] = raw_pdb
                st.session_state["current_pdb"] = target_info['pdb']
                st.success("Proteína baixada, limpa e preparada como PDBQT!")
            else:
                st.error("Erro na conversão para PDBQT.")
        else:
            st.error("Falha ao se conectar com o RCSB PDB.")

if "protein_data" in st.session_state:
    with st.expander("Inspeção Estrutural 3D (Proteína + Caixa)"):
        show_co_crystal = st.checkbox("Mostrar Ligante Co-cristalizado (Redocking)", value=True, key="co_crystal_1")
        view = py3Dmol.view(width=800, height=400)
        view.addModel(st.session_state["protein_data"], 'pdbqt')
        view.setStyle({'model': 0}, {'cartoon': {'color': 'spectrum'}, 'stick': {'radius': 0.1}})
        
        if show_co_crystal and "raw_pdb" in st.session_state:
            view.addModel(st.session_state["raw_pdb"], 'pdb')
            view.setStyle({'model': 1}, {}) # Esconder cadeias do raw_pdb
            view.setStyle({'model': 1, 'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2}})

        view.addBox({
            'center': {'x': center_x, 'y': center_y, 'z': center_z},
            'dimensions': {'w': size_x, 'h': size_y, 'd': size_z},
            'color': 'red', 'opacity': 0.5
        })
        view.zoomTo()
        showmol(view, height=400, width=800)

# --- Fase 2: Ligantes ---
st.header("2. Triagem de Ligantes (Diabetes)")
st.write("Insira os *SMILES* dos candidatos a inibidores que você quer testar contra este alvo:")
smiles_input = st.text_area("Exemplo:\nCC1=C(C=C(C=C1)C(=O)NC2=CC=C(C=C2)C3=NN(C(=N3)C)C)C", height=100)

ligands_to_dock = []
if smiles_input:
    for i, line in enumerate(smiles_input.strip().split('\n')):
        parts = line.strip().split()
        if parts:
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else f"Candidato_DB_{i+1}"
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mol = Chem.AddHs(mol)
                ligands_to_dock.append((name, mol))

if ligands_to_dock:
    st.info(f"{len(ligands_to_dock)} ligante(s) aceito(s) topologicamente.")

# --- Fase 3: Docking Pipeline ---
st.header("3. Executar Docking (AutoDock Vina)")
if st.button("Iniciar Simulação (Diabetes Pipeline)", type="primary"):
    if not ligands_to_dock:
        st.warning("Adicione pelo menos um ligante SMILES válido para prosseguir.")
    else:
        results = []
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        # Temp file for receptor
        t_rec = tempfile.NamedTemporaryFile(delete=False, suffix=".pdbqt")
        t_rec.write(st.session_state["protein_data"].encode('utf-8'))
        t_rec.close()
        
        v = Vina(sf_name='vina')
        v.set_receptor(t_rec.name)
        v.compute_vina_maps(center=[center_x, center_y, center_z], box_size=[size_x, size_y, size_z])
        
        total = len(ligands_to_dock)
        for idx, (lname, lmol) in enumerate(ligands_to_dock):
            status_text.text(f"Otimizando 3D e docando: {lname} ({idx+1}/{total})")
            
            # Etapa RDKit/Meeko
            AllChem.EmbedMolecule(lmol, AllChem.ETKDGv3())
            try:
                AllChem.UFFOptimizeMolecule(lmol)
            except: pass
            
            meeko_prep = MoleculePreparation()
            meeko_prep.prepare(lmol)
            lig_pdbqt = meeko_prep.write_pdbqt_string()
            
            # Etapa Vina
            if lig_pdbqt:
                v.set_ligand_from_string(lig_pdbqt)
                v.dock(exhaustiveness=exhaustiveness, n_poses=5)
                energies = v.energies(n_poses=1)
                best_affinity = energies[0][0] if len(energies)>0 else 0.0
                
                t_out = tempfile.NamedTemporaryFile(delete=False, suffix=".pdbqt")
                v.write_poses(t_out.name, n_poses=1, overwrite=True)
                with open(t_out.name, "r") as f:
                    out_pdbqt = f.read()
                os.unlink(t_out.name)
                
                results.append({"Ligand": lname, "Affinity (kcal/mol)": best_affinity, "Pose PDBQT": out_pdbqt})
            
            progress_bar.progress((idx+1)/total)
        
        os.unlink(t_rec.name)
        status_text.success("Pipeline da Diabetes finalizado!")
        
        df_res = pd.DataFrame(results)
        st.dataframe(df_res[["Ligand", "Affinity (kcal/mol)"]].sort_values("Affinity (kcal/mol)"))
        
        if not df_res.empty:
            st.write("### Melhor Complexo (Virtual Screening)")
            show_co_crystal_res = st.checkbox("Sobrepor Ligante Original (Redocking)", value=True, key="co_crystal_2")
            
            view_res = py3Dmol.view(width=800, height=500)
            view_res.addModel(st.session_state["protein_data"], 'pdbqt')
            view_res.setStyle({'model': 0}, {'cartoon': {'color': 'white'}})
            
            best_lig = df_res.iloc[0]["Pose PDBQT"]
            view_res.addModel(best_lig, 'pdbqt')
            view_res.setStyle({'model': 1}, {'stick': {'colorscheme': 'magentaCarbon', 'radius': 0.2}})
            
            if show_co_crystal_res and "raw_pdb" in st.session_state:
                view_res.addModel(st.session_state["raw_pdb"], 'pdb')
                view_res.setStyle({'model': 2}, {}) 
                view_res.setStyle({'model': 2, 'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2}})

            view_res.zoomTo()
            showmol(view_res, height=500, width=800)
