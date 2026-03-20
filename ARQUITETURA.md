# Arquitetura da Plataforma de Docking Molecular

Este documento estabelece a base de arquitetura de software, as bibliotecas recomendadas e o fluxo de dados projetado para a nova plataforma de docking molecular. Nosso objetivo é proporcionar um sistema seguro, modular, escalável e de fácil utilização (web-based) para pesquisadores e cientistas.

---

## 1. Bibliotecas e Tecnologias Sugeridas

A stack de tecnologia focará puramente no ecossistema de Python, com bibliotecas essenciais para web, bioinformática, quiminformática e computação científica de simulação.

### 1.1 UI e Framework Web
- **Streamlit**: Framework central que servirá o frontend do projeto. Ele permite escrever em Python scripts completos e rapidamente criar componentes reativos (inputs, barras laterais, e botões). Devido ao uso acadêmico de ML em python, é hoje preferido por dispensar conhecimento complexo de JS/React/Vue por parte dos cientistas.
- **stmol** / **py3Dmol**: Módulos vitais do Streamlit / Python para visualização ágil orientada a navegadores. Renderizarão a proteína macromolecular junto de seu grid box interativo e as exibições pós-docking do ligante na sua cavidade.
- **Plotly** / **Matplotlib**: Usados para a geração de gráficos de histogramas e estatísticas comparativas durante estudos de *Virtual Screening*.

### 1.2 Manipulação Química e Preparação Múltipla
- **RDKit**: O canivete suíço em quiminformática. Utilizado amplamente no backend do projeto para ingerir strings 2D (SMILES), arquivos SDF, realizar higienização química, gerar geometrias 3D iniciais com otimização (funções de campo de força como UFF ou MMFF94) antes de as passar adiante.
- **OpenBabel / Pybel**: Wrapper de biblioteca C++. Complemento poderoso frequentemente utilizado para lidar com a conversão massiva de formatos.
- **Biopython**: Manipular o arquivo de receptor `.pdb`, ajudando em lógicas de extração de sequências, limpezas, entre outras inspeções de resíduos necessárias antes da inserção no espaço virtual da simulação.

### 1.3 Motores de Preparação Avançada a Docking (Conversão Mágica PDBQT)
- **Meeko**: Moderno pacote python apoiado no próprio centro formador do AutoDock. Sua principal função na nossa pipeline é a criação do `ligand.pdbqt` direto de instâncias do RDKit sem depender do arriscado e obsoleto binário MGLTools.
- **PDB2PQR / ADFRsuite (reduce)**: Ajudantes potentes no momento da preparação da proteína para adicionar os hidrogênios do modo mais biologicamente favorável perante pH de simulações. Prepara-se em `receptor.pdb` passando para *pdbqt*.

### 1.4 Core de Simulação de Avaliação 
- **AutoDock Vina** (ou seu wrapper moderno em Python, a biblioteca `vina`): Motor consolidado globalmente. Trás robustez nos resultados em pouquíssimos segundos pelo uso do processamento multicore. Ele receberá a caixa dimensional, o receptor rígido e o ligante perfeitamente flexível.
- **Pandas**: Para agregar tabelas de `Score (kcal/mol)` por afinidades de energia livre, manipulação de RMSD correspondente as poses geradas do output `vina`.

---

## 2. Estrutura de Pastas e Diretórios (Modular)

Essa arquitetura tenta desacoplar a "apresentação" (`ui`) da "inteligência e validação" (`core`), fornecendo fáceis caminhos de manutenção.

```text
plataforma_docking/
│
├── app.py                      # Arquivo principal (Ponto de entrada: `streamlit run app.py`)
├── requirements.txt            # Dependências Python (pip)
├── README.md                   # Instrução de inicialização do repositório
├── ARQUITETURA.md              # Este documento guia de estrutura
│
├── core/
│   ├── __init__.py             # Módulo
│   ├── protein_prep.py         # Recebe .PDB, realiza limpeza, add/remove solvente, calcula carga e lança PDBQT
│   ├── ligand_prep.py          # Otimização 3D pelo RDKit e transformação para PDBQT flexível usando Meeko
│   ├── docking_engine.py       # Interação com a API Python do AutoDock Vina, definindo grid_center e volume
│   └── post_processing.py      # Extração das poses, parsing do log de score em dict/Dataframe
│
├── ui/
│   ├── __init__.py
│   ├── sidebar.py              # Concentra lógica e widgets da side-bar do Streamlit
│   ├── render_3d.py            # Trabalha especificamente o stmol / py3Dmol
│   └── tables.py               # Visualização dos DataFrames formatados no frontend
│
├── utils/
│   ├── filesystem_utils.py     # Geradores de ID seguros, pastas temporárias isoladas por sessão do usuário
│   └── logger.py               # Logs internos do processo do Backend
│
└── workspace/                  # Diretório transiente isolado, excluído do git via '.gitignore'
    ├── temp/                   # Local seguro global pra buffer
    └── sessions/               # Isolamento por sessão (Evita colisão de dados entre Múltiplos Usuários simultâneos)
        └── session_abc123/
            ├── inputs/         # arquivos virgens enviados pelo usuario e também os pre-processados .pdbqt
            └── outputs/        # Posições de resultado de docagem e seus logs (.csv e _out.pdbqt)
```

---

## 3. Fluxo Funcional de Vida (Life Cycle)

Este é o desenho passo a passo do trajeto de tela e de processamento de um pedido simulado da entrada à saída. Cada passo é reativo sob os pilares do ecossistema assíncrono do front para as chamadas sincronizadas subjacentes.

### Passo 1: Ingestão / Entrada (Interface Inicial)
1. Ao acessar a URL da plataforma, é gerado um ID de sessão única virtual (`session_abc123`) com seus diretórios privados.
2. Na interface (`sidebar.py` ou central), o usuário adiciona:
   - **Receptor:** O arquivo PDB ou PDBQT direto. 
   - **Ligantes:** Insere *SMILES*, `* .sdf` ou `*.mol2`. 
3. Os dados brutos injetados via `Streamlit Upload Widget` são recebidos por bytearrays e salvos pela camada de `filesystem_utils` no diretório `/workspace/sessions/<id>/inputs/`.

### Passo 2: Pré-processamento e Harmonização de Estruturas (O pipeline invisível)
1. **Tratativa Ligante** -> `core/ligand_prep.py` é acordado. Realiza uma checagem geométrica do 3D usando RDKit (congelando otimizações MMFF94) e gera o arquivo `ligante.pdbqt` pronto com o Meeko, reconhecendo toda rotatividade de torção das ligações.
2. **Tratativa Proteína** -> `core/protein_prep.py` tira resíduos de água inúteis do PDB caso precise, arrasta hidrogênios e cargas essenciais, gerando e salvando a saída `receptor.pdbqt`.

### Passo 3: Parametrização da Caixa-Alvo (Grid Box Targeting)
1. Para a simulação, é desenhado via componente Web (`ui/render_3d.py`) a proteína no espaço 3D.
2. O Usuário é instigado a ajustar cursores de `Centro [X,Y,Z]` e o tamanho da Caixa Virtual `Tamanhos [sX, sY, sZ]` que delimita as fronteiras moleculares ativas até englobar esteticamente a bolsa ativa do receptor na Web. Esses metadados numéricos do volume de simulação são retidos na memória (via `st.session_state`).

### Passo 4: Execução do Docking Físico
1. A ação é disparada pelo evento de Botão "Rodar Docking".
2. No Backend, `core/docking_engine.py` isola a máquina física. Com os dados extraídos: (*coordenadas grid box*, caminho do *receptor.pdbqt*, e *ligante.pdbqt*), ele instancia o motor **AutoDock Vina** (ou biblioteca Python).  
3. O software dispara o cálculo, o frontend desenha um visual de carregando (spinner). Os tensores varrem exaustivamente as congregações espaciais e a função de custo energética devolve uma lista de afinações preditivas ordenada, chamada "Poses".

### Passo 5: Exposição dos Resultados Analíticos e Download
1. As simulações terminadas formam novos `.pdbqt` em `/workspace/sessions/<id>/outputs/`.
2. O sistema de análise (`core/post_processing.py`) divide os scores de "Kcal/Mol".  
3. **Visão Tabular**: A interface pinta as métricas via tabela interativa `ui/tables.py` contendo também a similaridade estrutural RMSD.
4. **Visão Bioquímica (3D)**: Ao selecionar da tabela a "Pose 1" por exemplo, a interface 3D sobrepõe a macóromolecula junto com a representação molecular dockada daquela pose em tempo real, permitindo aos usuários avaliarem com minuciência os contatos hidrofóbicos/hidrofílicos formados.
5. **Encerramento de ciclo**: Um botão flutuante propõe ao usuário "Baixar Pacote de Dados" (.zip em stream), trazendo o PDB finalizado e as tabelas com os scores exportadas perfeitamente. Nas chamadas de encerramento do script, as tabelas temporárias em `workspace/sessions/` correspondendo à quele usuário são destruídas de modo ecológico pelo HD e RAM.
