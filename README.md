# Estimação dos paramêtros da distribuição NadaraJah-Naghighi
#### Objetiva-se avaliar numericamente estimadores da distribuição Nadarajah-Naghighi por meio de simulações de Monte Carlo. 

### 1. Introdução
O presente estudo aborda principalmente distribuições da área de sobrevivência. De forma geral, a literatura define como o tempo até a ocorrência de um evento (descrito por uma variável resposta) denota-se por tempo de sobrevivência ou tempo de falha e o conjunto de observações são os dados de sobrevivência.
Estes possuem suporte apenas nos reais positivos, e portanto, a distribuição normal não se caracteriza como adequada para descrevê-los.

A distribuição exponencial é amplamente empregada para problemas de análise de sobrevivência (área da saúde) ou confiabilidade (pesquisa industrial). Por exemplo, possui extensa aplicação na modelagem do tempo de vida de diversos produtos e materiais.	O foco deste relatório é voltada para a distribuição de Nadarajah-Haghighi (NH), proposta como uma generalização da distribuição exponencial. Sua função densidade de probabilidade é dada por:

$$ g(x)=\alpha \lambda(1+\lambda x)^{\alpha-1} \exp \left[ 1-\left( 1+\lambda x\right) ^\lambda\right], \quad \lambda>0,\quad \alpha>0, $$

em que $\lambda$ é o parâmetro de escala e $\alpha$ representa o parâmetro de forma. Quando $\alpha=1$  temos o caso particular da distribuição NH  equivalente à distribuição exponencial.

<!--- ![Fig1](https://raw.githubusercontent.com/carlosmanchini/parameter_estimation_nadarajah/fig1.png) --->
![Fig1](./main/fig1.png)

$$ G(x)=1-\exp \left[ 1-(1+\lambda x)^\alpha \right]. $$

A função de taxa de risco ou taxa de falha é uma quantidade que caracteriza o comportamento da sobrevivência até um determinando tempo. Em geral, a taxa de falha pode possuir um crescimento ou decrescimento não-monótono. 	Analisando-a pode-se definir o modelo probabilístico mais adequado para modelar o tempo de sobrevivência.	 Destaca-se que a distribuição NH é capaz de modelar dados que possuem taxas de falha não-monótonas independente da variação da densidade da distribuição. Este é um importante diferencial sob distribuições usuais em análise de sobrevivência como gama, Weibull e exponencial exponencializada que apenas acomodam taxas de falhas crescentes/decrescentes quando sua respectiva densidade estão diminuindo monotonicamente.	 
	A taxa de risco para a distribuição NH é expressa por: 
  
  $$ 	h(x) = \alpha \lambda (1+\lambda x)^{\alpha-1}. $$
  
  Generalizando uma distribuição é possível obter uma função de falha mais flexível, podendo ser constante, decrescente, crescente, em forma da banheira e banheira invertida, a depender dos valores dos parâmetros.
  
  O presente trabalho visa uma avaliação numérica dos EMV (Estimadores de Máxima Verossimilhança) $\hat{\alpha}$ e $\hat{\lambda}$ da distribuição NH, via simulação Monte Carlo. O comportamento e desempenho dos EMV será comparado sob amostras de tamanho finito e sob a inclusão do vetor escore com primeiras derivadas analíticas no processo de otimização de máxima verossimilhança, os quais serão denotados por $\hat{\alpha}_s$ e $\hat{\lambda}_s$.

### 2. Estimação Pontual
