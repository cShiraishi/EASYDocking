# EASYDocking 🧬

Esta aplicação permite realizar docking molecular utilizando **Streamlit** e **AutoDock Vina**.

## Requisitos de Sistema

-   Anaconda ou Miniconda instalado.
-   Python 3.8 ou superior.

## Instalação

1.  **Instale as dependências Python** (exceto Vina):

    ```bash
    pip install -r requirements.txt
    ```

2.  **Instale o AutoDock Vina via Conda** (Necessário pois a compilação via pip pode falhar no Mac):

    ```bash
    conda install -c conda-forge vina
    ```

    *Nota: Se o comando `conda` não estiver disponível, tente instalar as bibliotecas de desenvolvimento do Boost e usar `pip install vina`, ou use um gerenciador de pacotes como Homebrew (`brew install vina`) e certifique-se que o executável está no PATH (embora a aplicação use a biblioteca Python `vina`, que requer instalação via pip ou conda).*

    **Recomendação Forte:** Use o comando `conda install -c conda-forge vina` dentro do seu ambiente.

## Execução

Para iniciar a plataforma:

```bash
streamlit run app.py
```

## Como Usar

1.  **Proteína**: Carregue um arquivo `.pdbqt` (do AutoDockTools/MGLTools) ou `.pdb`.
2.  **Ligantes**: Cole SMILES ou carregue um `.sdf`.
3.  **Grid**: Ajuste o centro e tamanho da caixa visualmente até cobrir o sítio ativo.
4.  **Docking**: Clique em "Iniciar Docking".
5.  **Resultados**: Visualize a melhor pose e baixe a tabela de afinidades.
