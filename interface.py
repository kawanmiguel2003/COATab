import pandas as pd
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import mplcursors
from scipy.stats import shapiro
from scipy.stats import entropy

# Configuração da página para largura total
st.set_page_config(layout="wide")

# Adicionando um título (cabeçalho) abaixo da imagem
st.markdown("<h1 style='text-align: center;'>COA - Tab</h1>", unsafe_allow_html=True)
st.markdown("<h3 style='text-align: center;'>Análise de Capacidade do Processo</h3>", unsafe_allow_html=True)

# Títulos e descrições para as abas
aba1, aba2, aba3 = st.tabs(["Entrada de Dados", "Resumo das Métricas", "Gráficos"])

with aba1:
    st.markdown("**Faça o upload de um arquivo CSV ou Excel para iniciar a análise.**")
    
    uploaded_file = st.file_uploader("Escolha um arquivo para análise:", type=["csv", "xlsx"])
    
    if uploaded_file:
        # Exibir o arquivo carregado
        if uploaded_file.name.endswith('.csv'):
            df = pd.read_csv(uploaded_file)
        else:
            df = pd.read_excel(uploaded_file)
        
        

        column = st.selectbox("Selecione a coluna para análise:", df.columns)

        if column:
            # Tentar converter a coluna para numérico
            df[column] = pd.to_numeric(df[column], errors='coerce')

            # Tratar valores ausentes
            if df[column].isna().sum() > 0:
                st.warning("Foram encontrados valores não numéricos ou ausentes. Eles serão ignorados.")
                df = df.dropna(subset=[column])

            # Inputs para Limites de Especificação com descrições
            st.markdown("**Defina os Limites de Especificação:**")
            USL = st.number_input("Limite Superior de Especificação (USL):", value=3.0, step=0.1, help="Valor máximo permitido para o processo.")
            LSL = st.number_input("Limite Inferior de Especificação (LSL):", value=1.0, step=0.1, help="Valor mínimo permitido para o processo.")

            st.dataframe(df)
            
            # Calcular as métricas
            metric = df[column].dropna()  # Remover NaN
            mean_metric = np.mean(metric)
            std_metric = np.std(metric, ddof=1)

            # Calcular Cp, Cpk, Pp, Ppk, Cpm, CpG, VPC, DPMO, PPI, Skewness e Kurtosis
            Cp = (USL - LSL) / (6 * std_metric)
            Cpk = min((USL - mean_metric) / (3 * std_metric), (mean_metric - LSL) / (3 * std_metric))
            Pp = (USL - LSL) / (6 * np.std(df[column], ddof=1))
            Ppk = min((USL - np.mean(df[column])) / (3 * np.std(df[column], ddof=1)), 
                      (np.mean(df[column]) - LSL) / (3 * np.std(df[column], ddof=1)))
            Cpm = (USL - LSL) / (6 * np.sqrt(std_metric**2 + (mean_metric - (USL + LSL) / 2)**2))
            CpG = (USL - LSL) / (6 * std_metric)
            VPC = (6 * std_metric) / (USL - LSL) * 100
            sigma_process = max((USL - mean_metric) / std_metric, (mean_metric - LSL) / std_metric)
            DPMO = (metric.isna().sum() / len(metric)) * 1_000_000
            PPI = (USL - LSL) / (6 * std_metric)
            skewness = metric.skew()
            kurtosis = metric.kurtosis()
            # Coeficiente de Variação (CV%)
            cv_percent = (std_metric / mean_metric) * 100


            # Entropia dos dados
            data_entropy = entropy(np.histogram(metric, bins=20)[0])

            # Intervalo de Confiança para a Média
            confidence_interval = stats.t.interval(0.95, len(metric)-1, loc=mean_metric, scale=stats.sem(metric))

            # Teste de Normalidade (Shapiro-Wilk)
            stat, p_value = shapiro(metric)

            # Colocar as métricas na aba 2 apenas quando o arquivo for carregado
            with aba2:
                col1, col2 = st.columns(2)

                with col1:
                    # Exibir os resultados
                    st.write(f"**Média:** {mean_metric:.4f}")
                    st.write(f"**Desvio Padrão:** {std_metric:.4f}")
                    st.write(f"**Cp:** {Cp:.4f}")
                    st.write(f"**Cpk:** {Cpk:.4f}")
                    st.write(f"**Pp:** {Pp:.4f}")
                    st.write(f"**Ppk:** {Ppk:.4f}")
                    st.write(f"**Cpm:** {Cpm:.4f}")
                    st.write(f"**CV%:** {cv_percent:.4f}")
                    st.write(f"**Entropia dos dados:** {data_entropy:.4f}")
                    st.write(f"**Intervalo de Confiança para a Média:** ({confidence_interval[0]:.4f}, {confidence_interval[1]:.4f})")
                with col2:
                    st.write(f"**CpG:** {CpG:.4f}")
                    st.write(f"**VPC:** {VPC:.2f}%")
                    st.write(f"**Sigma do Processo:** {sigma_process:.4f}")
                    st.write(f"**DPMO:** {DPMO:.0f}")
                    st.write(f"**PPI (Índice de Performance):** {PPI:.4f}")
                    st.write(f"**Skewness:** {skewness:.4f}")
                    st.write(f"**Kurtosis:** {kurtosis:.4f}")
                    st.write(f"**Shapiro-Wilk Test Stat:** {stat:.4f}")
                    st.write(f"**Shapiro-Wilk p-value:** {p_value:.4f}")

with aba3:
    # Aqui, os gráficos de análise continuam, mas apenas se o arquivo for carregado
    if uploaded_file:
        # Gráfico Q-Q (Quantile-Quantile)
        fig, ax = plt.subplots(figsize=(15, 6))
        stats.probplot(metric, dist="norm", plot=ax)
        ax.set_title(f"Gráfico Q-Q da Coluna {column}")
        st.pyplot(fig)

        # Adicionar separador após o gráfico
        st.markdown("<hr>", unsafe_allow_html=True)
        
        # Calcular a amplitude móvel
        window_size = st.slider("Escolha o tamanho da janela para o gráfico de amplitude móvel:", 2, 50, 10)
        rolling_amplitude = metric.rolling(window=window_size).apply(lambda x: x.max() - x.min())
        
        # Cálculos da média e limites de controle
        media_amplitude_movel = np.mean(rolling_amplitude.dropna())
        limite_superior_controle = media_amplitude_movel + 3 * std_metric
        limite_inferior_controle = media_amplitude_movel - 3 * std_metric

        # Gráfico de Amplitude Móvel
        fig, ax = plt.subplots(figsize=(15, 6))
        ax.plot(rolling_amplitude, color='blue', linestyle='-', marker='o', markersize=5, label='Amplitude Móvel')
        ax.axhline(y=media_amplitude_movel, color='green', linestyle='--', linewidth=1.5, label=f'Média Amp.Móvel = {media_amplitude_movel:.3f}')
        ax.axhline(y=limite_superior_controle, color='red', linestyle='--', linewidth=1.5, label=f'Limite Superior = {limite_superior_controle:.3f}')
        ax.axhline(y=limite_inferior_controle, color='red', linestyle='--', linewidth=1.5, label=f'Limite Inferior = {limite_inferior_controle:.3f}')
        ax.set_title('Gráfico de Amplitude Móvel', fontsize=12)
        ax.set_xlabel('Índice', fontsize=10)
        ax.set_ylabel('Amplitude', fontsize=10)
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        mplcursors.cursor()
        st.pyplot(fig)
        
        # Adicionar separador após o gráfico
        st.markdown("<hr>", unsafe_allow_html=True)

        # Gráfico de Valores Individuais
        fig, ax = plt.subplots(figsize=(15, 6))
        ax.scatter(range(len(metric)), metric, color='blue', edgecolor='black', label='Valores Individuais', s=50)
        ax.axhline(y=mean_metric, color='green', linestyle='--', linewidth=2, label=f'Média = {mean_metric:.3f}')
        ax.axhline(y=USL, color='blue', linestyle='--', linewidth=2, label=f'USL = {USL:.3f}')
        ax.axhline(y=LSL, color='red', linestyle='--', linewidth=2, label=f'LSL = {LSL:.3f}')
        ax.set_title(f'Valores Individuais para {column}', fontsize=12)
        ax.set_xlabel('Índice', fontsize=10)
        ax.set_ylabel(f'{column}', fontsize=10)
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        ax.grid(True, axis='y')
        mplcursors.cursor()
        st.pyplot(fig)
        
        # Adicionar separador após o gráfico
        st.markdown("<hr>", unsafe_allow_html=True)

        # Gráfico de Boxplot com Limites de Especificação
        fig, ax = plt.subplots(figsize=(15, 6))
        sns.boxplot(data=metric, ax=ax, color="lightblue", width=0.5)
        ax.axhline(y=USL, color='red', linestyle='--', linewidth=2, label=f'USL ({USL})')
        ax.axhline(y=LSL, color='green', linestyle='--', linewidth=2, label=f'LSL ({LSL})')
        ax.set_title(f"Boxplot de {column} com Limites de Especificação", fontsize=12)
        ax.set_xlabel(f"{column}", fontsize=10)
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        mplcursors.cursor()
        st.pyplot(fig)

        # Controle deslizante para ajustar o número de bins
        bins = st.slider(
            "Escolha o número de bins para o histograma:", 5, 100, 30, step=1
        )
        # Gráfico de Capabilidade do Processo
        fig, ax = plt.subplots(figsize=(15, 6))
        sns.histplot(metric, kde=True, bins=bins, ax=ax, color='blue')
        ax.axvline(mean_metric, color='green', linestyle='--', linewidth=2, label=f'Média = {mean_metric:.2f}')
        ax.axvline(USL, color='red', linestyle='--', linewidth=2, label=f'USL = {USL:.2f}')
        ax.axvline(LSL, color='red', linestyle='--', linewidth=2, label=f'LSL = {LSL:.2f}')
        ax.set_title(f"Histograma de {column} com Limites de Especificação", fontsize=12)
        ax.set_xlabel(f'{column}', fontsize=10)
        ax.set_ylabel("Frequência", fontsize=10)
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        st.pyplot(fig)

        
        # Adicionar separador após o gráfico
        st.markdown("<hr>", unsafe_allow_html=True)

        # Gráfico X-barra e R
        fig, axs = plt.subplots(2, 1, figsize=(20, 12), sharex=True)
        axs[0].plot(metric, marker='o', color='blue', label='Média dos Subgrupos')
        axs[0].axhline(y=mean_metric, color='green', linestyle='--', linewidth=2, label=f'Média ({mean_metric:.3f})')
        axs[0].axhline(y=mean_metric + 3*std_metric, color='red', linestyle='--', linewidth=2, label=f'UCL ({mean_metric + 3*std_metric:.3f})')
        axs[0].axhline(y=mean_metric - 3*std_metric, color='red', linestyle='--', linewidth=2, label=f'LCL ({mean_metric - 3*std_metric:.3f})')
        axs[0].set_title("Gráfico X-barra", fontsize=12)
        axs[0].set_ylabel(column, fontsize=10)
        axs[0].legend()
        mplcursors.cursor()

        r_values = metric.rolling(window=2).apply(lambda x: max(x) - min(x)).dropna()
        axs[1].plot(r_values, marker='o', color='blue', label='Amplitude dos Subgrupos')
        axs[1].axhline(y=r_values.mean(), color='green', linestyle='--', linewidth=2, label=f'Média ({r_values.mean():.3f})')
        axs[1].axhline(y=r_values.mean() + 3*std_metric, color='red', linestyle='--', linewidth=2, label=f'UCL ({r_values.mean() + 3*std_metric:.3f})')
        axs[1].axhline(y=r_values.mean() - 3*std_metric, color='red', linestyle='--', linewidth=2, label=f'LCL ({r_values.mean() - 3*std_metric:.3f})')
        axs[1].set_title("Gráfico R (Amplitude)", fontsize=12)
        axs[1].set_xlabel("Subgrupo", fontsize=10)
        axs[1].set_ylabel("Amplitude", fontsize=10)
        axs[1].legend()
        mplcursors.cursor()
        st.pyplot(fig)

