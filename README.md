# Un nuevo método de Monte Carlo para determinar la distribución nula empírica bajo la hipótesis de no asociación para GWAS

## Datos reales

Para este estudio se utilizó un conjunto de datos de la raza bovina de carne **Rubia Gallega**, proporcionado por la **Asociación Nacional de Criadores de Ganado Vacuno Selecto de Raza Rubia Gallega (ACRUGA)**. El conjunto de datos genealógico incluye un total de **541.492 registros individuo–padre–madre**.

El conjunto de datos genotípicos consta de **5.627 individuos**, genotipados mediante los arrays **Axiom_BovMDv2** (688 individuos) y **Axiom_BovMDv3** (4.939 individuos). Ambos arrays fueron procesados con la plataforma **Axiom Bovine (ThermoFisher Scientific)** con el fin de unificar mapas y genotipos. Tras el control de calidad, se retuvieron **39.194 SNPs** para el análisis.

El rasgo analizado fue el **peso de la canal en frío (Cold Carcass Weight, CCW)**. El conjunto fenotípico incluye **99.227 registros** procedentes de **41 mataderos**, correspondientes a animales con edades comprendidas entre **200 y 400 días**. El CCW presentó una media fenotípica de **224,99 kg** y una desviación estándar de **41,31 kg**.

Los animales fueron criados en **58 explotaciones**, bajo **1.860 condiciones productivas temporales distintas**, definidas mediante una clasificación año–estación.

---

## Modelo estadístico

La predicción de valores genéticos se realizó mediante un modelo **ssGBLUP (single-step Genomic Best Linear Unbiased Prediction)**, que puede expresarse en notación matricial como:

```
y = d + Xb + Tr + Wu + e
```

donde:

* **y** es el vector de fenotipos del rasgo CCW
* **d** es el vector de edad en días, incluido como covariable
* **b** es el vector de efectos fijos (sexo, matadero y explotación)
* **r** es el vector de efectos aleatorios año–estación
* **u** es el vector de valores genéticos aditivos
* **e** es el vector de efectos residuales
* **X**, **T** y **W** son las matrices de incidencia correspondientes

Los efectos aleatorios se distribuyeron como:

* r ~ N(0, Iσ²_R)
* u ~ N(0, Hσ²_A)
* e ~ N(0, Iσ²_E)

La matriz **H** corresponde a la matriz de relaciones genéticas single-step, que combina la información genealógica (**A**) y genómica (**G**, Van Raden, 2008), siguiendo a Christensen et al. (2009) y Aguilar et al. (2010). Su inversa se calculó como:

```
H⁻¹ = A⁻¹ +
       | 0              0 |
       | 0   G⁻¹ − A⁻¹₂₂ |
```

Las varianzas se fijaron en:

* Varianza año–estación (σ²_R): **195,9 kg²**
* Varianza genética aditiva (σ²_A): **65,3 kg²**
* Varianza residual (σ²_E): **522,4 kg²**

---

## Análisis estadísticos

### ssGBLUP y GWAS

Los valores genéticos fueron predichos utilizando el software **BLUPF90+** (Lourenco et al., 2022). Los efectos de los SNP (â) se estimaron con **POSTGSF90** (Misztal et al., 2014).

La varianza asociada al marcador *k* se calculó como:

```
Var(â_k) = 2 p_k (1 − p_k) â_k²
```

donde *p_k* es la frecuencia alélica del SNP *k*.

---

## Distribución estocástica de asociación

Una propiedad conocida del modelo **SNP-BLUP** es su equivalencia con **GBLUP** (Wang et al., 2012). Esta propiedad se utilizó para definir una distribución estocástica de asociación de SNPs mediante un algoritmo computacionalmente viable.

El efecto de sustitución alélica del marcador *k* (a_k) se muestreó de una distribución normal estándar:

```
a_k ~ N(0, 1)
```

A partir de estos efectos, el valor genético del individuo *i* se calculó como:

```
BV_i = Σ_k M_(i,k) · a_k
```

donde *M_(i,k)* es el valor del SNP *k* en el individuo *i* de la matriz de dosis alélicas centradas **M**.

En notación matricial:

```
BV = M ∘ (1_a · a')
```

donde ∘ representa el producto de Hadamard.

A partir de los valores genéticos simulados se obtuvieron efectos SNP y varianzas de marcador mediante GWAS (Zhang et al., 2010). Este procedimiento se repitió **1.000 veces** para construir la distribución estocástica de varianzas de asociación.

---

## Comparación entre GWAS y distribución estocástica

La varianza obtenida mediante GWAS se comparó con la distribución estocástica de asociación. En cada réplica, la varianza asociada a los marcadores se reescaló a la varianza genética total explicada por SNPs.

Siguiendo a Sahana et al. (2010) y Li et al. (2021), las varianzas se agruparon en **ventanas de 25 SNPs**, tanto para los resultados de GWAS como para la distribución estocástica.

Las varianzas explicadas por cada ventana se compararon entre ambos enfoques. Para cada ventana se calculó un **valor p empírico**, definido como la probabilidad de que la varianza explicada por GWAS pertenezca a la distribución estocástica obtenida mediante ssGBLUP.

Las ventanas cuya varianza explicada fue inferior a la esperada estocásticamente fueron descartadas. Los valores p se corrigieron mediante **corrección de Bonferroni**, considerando el número de tests independientes realizados.

---

## Análisis del array de ADN

A partir de las posiciones físicas (bp) del mapa genómico del array de ADN, se calculó la **densidad de SNPs en fragmentos no solapados de 1 Mb**. Estos resultados se compararon con los obtenidos en el GWAS y con los valores p empíricos derivados del análisis de asociación.
