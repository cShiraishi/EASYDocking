import streamlit as st
import tempfile
import os
from langchain_community.document_loaders import PyPDFLoader
from langchain_text_splitters import RecursiveCharacterTextSplitter
from langchain_huggingface import HuggingFaceEmbeddings
from langchain_community.chat_models import ChatOllama
from langchain_community.vectorstores import Chroma
from langchain_classic.chains import create_retrieval_chain
from langchain_classic.chains.combine_documents import create_stuff_documents_chain
from langchain_core.prompts import ChatPromptTemplate

# Configuração da página
st.set_page_config(page_title="Assistente de Artigos - RAG", layout="wide")

st.title("Assistente RAG para Artigos Científicos 📚 (100% Local)")
st.markdown("Com esta ferramenta, você pode fazer upload de um artigo em PDF e conversar com ele. Todo o processamento e inferência ocorrem na sua máquina sem enviar dados para a internet, utilizando **Embeddings do HuggingFace** e **Llama 3 via Ollama**.")

# Variável global da sessão para armazenar o histórico e a base de vetores
if "messages" not in st.session_state:
    st.session_state.messages = []

if "vectorstore" not in st.session_state:
    st.session_state.vectorstore = None

st.sidebar.header("1. Upload de Documentos")
uploaded_file = st.sidebar.file_uploader("Selecione o artigo (PDF)", type="pdf")

if uploaded_file is not None:
    # Apenas processa se for um novo PDF para evitar recálculos constantes
    if st.sidebar.button("Processar Documento e Inicializar RAG"):
        st.session_state.messages = []
        
        # Salva PDF no diretório temporário
        temp_dir = tempfile.mkdtemp()
        path = os.path.join(temp_dir, uploaded_file.name)
        with open(path, "wb") as f:
            f.write(uploaded_file.getvalue())

        with st.spinner("Lendo documento... (1/3)"):
            loader = PyPDFLoader(path)
            docs = loader.load()

        with st.spinner("Criando Chunks de texto... (2/3)"):
            # Otimização de chunking 
            text_splitter = RecursiveCharacterTextSplitter(chunk_size=1000, chunk_overlap=200)
            splits = text_splitter.split_documents(docs)

        with st.spinner("Baixando/Carregando embeddings locais e vetorizando... (3/3 - Na 1ª vez pode demorar para baixar)"):
            # Embeddings locais (SentenceTransformers - grátis e offline)
            embeddings = HuggingFaceEmbeddings(model_name="all-MiniLM-L6-v2")
            vectorstore = Chroma.from_documents(documents=splits, embedding=embeddings)
            st.session_state.vectorstore = vectorstore
            
        st.sidebar.success("Documento carregado e vetorizado no ChromaDB!")

# Exibição do Histórico do Chat
for message in st.session_state.messages:
    with st.chat_message(message["role"]):
        st.markdown(message["content"])

# Só mostra o input de conversa quando tem documento
if st.session_state.vectorstore is not None:
    # Caixa principal para o novo prompt
    user_question = st.chat_input("2. Pergunte algo sobre os métodos apontados no artigo:")
    
    if user_question:
        # Adiciona no chat visualmente e na sessão
        st.session_state.messages.append({"role": "user", "content": user_question})
        with st.chat_message("user"):
            st.markdown(user_question)
            
        # Gera a resposta via LLM
        with st.chat_message("assistant"):
            with st.spinner("Llama 3 está pensando... (Via Ollama)"):
                try:
                    # Carregando retriever do vector database usando a busca k=3
                    retriever = st.session_state.vectorstore.as_retriever(search_kwargs={"k": 3})
                    
                    # Carregando Llama 3 via Ollama
                    llm = ChatOllama(model="llama3", temperature=0.1)
                    
                    # Define prompt para evitar "alucinação" (inventar respostas fora do texto)
                    system_prompt = (
                        "Você é um excelente assistente científico de suporte à leitura. "
                        "Sua principal responsabilidade é responder o usário com base no artigo lido."
                        "Utilize APENAS os contextos retirados do artigo listados abaixo para responder à pergunta."
                        "Se você não souber ou não encontrar a resposta nos trechos de contexto, "
                        "simplesmente diga 'Desculpe, com base no conteúdo lido, não localizei a resposta para isso.'\n\n"
                        "Contexto do Artigo:\n{context}"
                    )
                    
                    prompt = ChatPromptTemplate.from_messages([
                        ("system", system_prompt),
                        ("human", "{input}"),
                    ])
                    
                    question_answer_chain = create_stuff_documents_chain(llm, prompt)
                    rag_chain = create_retrieval_chain(retriever, question_answer_chain)
                    
                    # Invocação
                    response = rag_chain.invoke({"input": user_question})
                    
                    answer_text = response["answer"]
                    st.markdown(answer_text)
                    st.session_state.messages.append({"role": "assistant", "content": answer_text})
                    
                    # Exibe fontes como extra no final da resposta
                    with st.expander("Trechos Consultados do Artigo (Fontes de Contexto)"):
                        for i, doc in enumerate(response["context"]):
                            page = doc.metadata.get('page', 0) + 1
                            st.info(f"**Referência {i+1} (Página {page}):**\n\n _{doc.page_content}_")
                
                except Exception as e:
                    st.error(f"Erro ao conversar com Llama3/Ollama: {e}")
                    st.warning("Verifique se o comando de pull do Llama 3 no fundo já finalizou.")
else:
    st.info("← Para iniciar, faça upload do documento PDF na barra lateral e clique em Processar.")
