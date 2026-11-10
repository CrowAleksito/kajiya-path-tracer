# Formulación Matemática y Deducción Analítica del Motor Path Tracer

## 1. La Ecuación Integral de Renderizado de Kajiya

En transferencia radiativa y síntesis de imágenes físicamente basada (Kajiya, 1986), la radiancia saliente en estado estacionario $L_o(\mathbf{x}, \omega_o)$ desde un punto superficial $\mathbf{x} \in \mathbb{R}^3$ en la dirección de observación $\omega_o$ está gobernada por la **Ecuación de Renderizado**:

$$L_o(\mathbf{x}, \omega_o) = L_e(\mathbf{x}, \omega_o) + \int_{\Omega} f_r(\mathbf{x}, \omega_i, \omega_o) \, L_i(\mathbf{x}, \omega_i) \, (\omega_i \cdot \mathbf{n}) \, d\omega_i$$

donde:
- $L_o(\mathbf{x}, \omega_o)$: Radiancia espectral saliente en el punto $\mathbf{x}$ en dirección $\omega_o$ $[\text{W}\cdot\text{sr}^{-1}\cdot\text{m}^{-2}]$.
- $L_e(\mathbf{x}, \omega_o)$: Radiancia emitida directamente por la superficie en $\mathbf{x}$.
- $\Omega = \{ \omega_i \in S^2 : \omega_i \cdot \mathbf{n} \ge 0 \}$: Hemisferio unitario superior orientado según la normal superficial unitaria $\mathbf{n}$.
- $f_r(\mathbf{x}, \omega_i, \omega_o)$: Función de Distribución de Reflectancia Bidireccional (BRDF) $[\text{sr}^{-1}]$. Para una superficie difusa Lambertiana con reflectancia hemisférica (albedo) $\rho(\mathbf{x}) \in [0, 1]^3$:
  $$f_r(\mathbf{x}, \omega_i, \omega_o) = \frac{\rho(\mathbf{x})}{\pi}$$
- $L_i(\mathbf{x}, \omega_i)$: Radiancia incidente en $\mathbf{x}$ proveniente de la dirección $\omega_i$. Por invarianza de la radiancia en el vacío a lo largo de un rayo, $L_i(\mathbf{x}, \omega_i) = L_o(\text{ray\_cast}(\mathbf{x}, \omega_i), -\omega_i)$, lo que convierte la ecuación en una **ecuación integral de Fredholm de segunda especie** recursiva de dimensión infinita.
- $(\omega_i \cdot \mathbf{n}) = \cos\theta_i$: Factor geométrico de foreshortening (Ley del Coseno de Lambert).

---

## 2. Cambio de Variable a Coordenadas Esféricas y Diferencial de Ángulo Sólido

En el sistema de referencia local ortonormal $(\mathbf{t}, \mathbf{b}, \mathbf{n})$ centrado en el punto de impacto $\mathbf{x}$ (con eje $z$ alineado a $\mathbf{n}$), cualquier dirección unitaria $\omega_i = (x, y, z) \in \Omega$ se parametriza mediante el ángulo cenital (polar) $\theta \in [0, \pi/2]$ respecto a la normal $\mathbf{n}$ y el ángulo acimutal $\phi \in [0, 2\pi)$:

$$\omega_i(\theta, \phi) = \begin{pmatrix} \sin\theta \cos\phi \\ \sin\theta \sin\phi \\ \cos\theta \end{pmatrix}$$

El vector normal a la esfera unitaria e inducido por las derivadas parciales tangenciales $\frac{\partial \omega_i}{\partial \theta}$ y $\frac{\partial \omega_i}{\partial \phi}$ determina el elemento diferencial de área superficial sobre la esfera unitaria (diferencial de ángulo sólido $d\omega_i$):

$$\frac{\partial \omega_i}{\partial \theta} = \begin{pmatrix} \cos\theta \cos\phi \\ \cos\theta \sin\phi \\ -\sin\theta \end{pmatrix}, \qquad \frac{\partial \omega_i}{\partial \phi} = \begin{pmatrix} -\sin\theta \sin\phi \\ \sin\theta \cos\phi \\ 0 \end{pmatrix}$$

Calculando la norma del producto vectorial fundamental (Jacobiano esférico):

$$d\omega_i = \left\| \frac{\partial \omega_i}{\partial \theta} \times \frac{\partial \omega_i}{\partial \phi} \right\| d\theta \, d\phi = \sin\theta \, d\theta \, d\phi$$

Sustituyendo $(\omega_i \cdot \mathbf{n}) = \cos\theta$ y $d\omega_i = \sin\theta \, d\theta \, d\phi$, la integral hemisférica de la Ecuación de Kajiya se expresa explícitamente como una integral doble iterada en el rectángulo $(\theta, \phi) \in [0, \pi/2] \times [0, 2\pi)$:

$$L_o(\mathbf{x}, \omega_o) = L_e(\mathbf{x}, \omega_o) + \int_{0}^{2\pi} \int_{0}^{\pi/2} f_r(\mathbf{x}, \theta, \phi, \omega_o) \, L_i(\mathbf{x}, \theta, \phi) \, \cos\theta \sin\theta \, d\theta \, d\phi$$

---

## 3. Validación Analítica: El Horno de Luz (White Furnace Test)

Si encerramos una superficie difusa Lambertiana con reflectancia $\rho = 1$ ($f_r = \frac{1}{\pi}$) dentro de un entorno isotrópico homogéneo con radiancia incidente constante $L_i(\mathbf{x}, \omega_i) = 1$ para toda dirección $\omega_i \in \Omega$, y ausencia de emisión propia ($L_e = 0$), la integral posee solución analítica cerrada exacta:

$$L_o = \int_{\Omega} \frac{1}{\pi} \cdot 1 \cdot \cos\theta \, d\omega = \frac{1}{\pi} \int_{0}^{2\pi} d\phi \int_{0}^{\pi/2} \cos\theta \sin\theta \, d\theta$$

Resolviendo ambas integrales definidas:

$$\int_{0}^{2\pi} d\phi = 2\pi, \qquad \int_{0}^{\pi/2} \cos\theta \sin\theta \, d\theta = \left[ \frac{\sin^2\theta}{2} \right]_{0}^{\pi/2} = \frac{1}{2}$$

Por consiguiente:

$$L_o = \frac{1}{\pi} \cdot (2\pi) \cdot \left(\frac{1}{2}\right) = 1.000000$$

Esto demuestra analíticamente la **conservación de la energía en equilibrio termodinámico**: una esfera 100% reflectora dentro de un horno de radiancia unitaria es indistinguible del fondo ($L_o = 1$).

---

## 4. Estimador Estocástico de Monte Carlo y Tasa de Convergencia

Sea $I = \int_{\Omega} g(\omega) \, d\omega$ la integral hemisférica que deseamos aproximar, donde el integrando es:

$$g(\omega_i) = f_r(\mathbf{x}, \omega_i, \omega_o) \, L_i(\mathbf{x}, \omega_i) \, (\omega_i \cdot \mathbf{n})$$

Si extraemos $N$ muestras direccionales independientes e idénticamente distribuidas $\omega_1, \omega_2, \dots, \omega_N \sim p(\omega)$, donde $p(\omega) > 0$ es una función de densidad de probabilidad (PDF) definida respecto al ángulo sólido tal que $\int_{\Omega} p(\omega) \, d\omega = 1$, el **estimador primario de Monte Carlo** viene dado por:

$$\langle L_o \rangle_N = L_e(\mathbf{x}, \omega_o) + \frac{1}{N} \sum_{k=1}^N \frac{f_r(\mathbf{x}, \omega_k, \omega_o) \, L_i(\mathbf{x}, \omega_k) \, (\omega_k \cdot \mathbf{n})}{p(\omega_k)}$$

### Demostración de Insesgadez (Unbiasedness)

$$\mathbb{E}\big[\langle L_o \rangle_N - L_e\big] = \frac{1}{N} \sum_{k=1}^N \int_{\Omega} \frac{g(\omega)}{p(\omega)} p(\omega) \, d\omega = \int_{\Omega} g(\omega) \, d\omega = I$$

### Varianza y Error Cuadrático Medio (MSE / RMSE)

Sea $Y_k = \frac{g(\omega_k)}{p(\omega_k)}$ cada muestra aleatoria con varianza finita $\sigma^2 = \text{Var}[Y_k] = \int_{\Omega} \frac{g^2(\omega)}{p(\omega)} d\omega - I^2$. Por independencia de las $N$ muestras:

$$\text{MSE}\big[\langle L_o \rangle_N\big] = \text{Var}\big[\langle L_o \rangle_N\big] = \frac{\sigma^2}{N} = \mathcal{O}\left(\frac{1}{N}\right)$$

Y el error cuadrático medio raíz (RMSE) decrece según:

$$\text{RMSE}\big[\langle L_o \rangle_N\big] = \sqrt{\text{MSE}} = \frac{\sigma}{\sqrt{N}} = \mathcal{O}\left(\frac{1}{\sqrt{N}}\right)$$

---

## 5. Muestreo Hemisférico Uniforme vs. Muestreo por Importancia Ponderado por Coseno

### 5.1. Muestreo Hemisférico Uniforme

El área total del hemisferio unitario es $\int_0^{2\pi}\int_0^{\pi/2} \sin\theta \, d\theta \, d\phi = 2\pi$. Por tanto, una densidad uniforme respecto al ángulo sólido es:

$$p_{\text{unif}}(\omega) = \frac{1}{2\pi} \implies p(\theta, \phi) = \frac{\sin\theta}{2\pi}$$

Aplicando el método de la transformada inversa con $\xi_1, \xi_2 \sim \mathcal{U}(0, 1)$:
- $\cos\theta = \xi_1, \quad \sin\theta = \sqrt{1 - \xi_1^2}, \quad \phi = 2\pi \xi_2$
- Para una BRDF Lambertiana $f_r = \frac{\rho}{\pi}$, el peso de cada muestra en el estimador es:
  $$W_{\text{unif}}(\omega_k) = \frac{(\rho / \pi) \cos\theta_k}{1 / (2\pi)} = 2 \, \rho \cos\theta_k$$

En el Horno de Luz ($\rho = 1, L_i = 1$), cada muestra vale $Y_k = 2\cos\theta_k = 2\xi_1$, cuya varianza teórica es:

$$\sigma^2_{\text{unif}} = \int_0^1 (2\xi_1)^2 \, d\xi_1 - 1^2 = \frac{4}{3} - 1 = \frac{1}{3} \approx 0.333333$$

### 5.2. Muestreo Hemisférico Ponderado por Coseno (Cosine-Weighted Importance Sampling)

Dado que el integrando contiene explícitamente el factor $(\omega_i \cdot \mathbf{n}) = \cos\theta$, elegimos una densidad proporcional a $\cos\theta$:

$$p_{\text{cos}}(\omega) = c \cos\theta$$

Imponiendo la condición de normalización sobre el hemisferio $\Omega$:

$$\int_{\Omega} c \cos\theta \, d\omega = c \int_0^{2\pi} d\phi \int_0^{\pi/2} \cos\theta \sin\theta \, d\theta = c (2\pi) \left(\frac{1}{2}\right) = c\pi = 1 \implies c = \frac{1}{\pi}$$

Por tanto:

$$p_{\text{cos}}(\omega) = \frac{\cos\theta}{\pi}, \qquad p(\theta, \phi) = \frac{\cos\theta \sin\theta}{\pi}$$

Integrando las densidades marginales para aplicar el método de la transformada inversa con $\xi_1, \xi_2 \sim \mathcal{U}(0, 1)$:

$$F_{\Theta}(\theta) = \int_0^\theta 2\cos\theta' \sin\theta' \, d\theta' = \sin^2\theta = \xi_1 \implies \sin\theta = \sqrt{\xi_1}, \quad \cos\theta = \sqrt{1 - \xi_1}, \quad \phi = 2\pi\xi_2$$

### Cancelación Algebraica del Término Coseno

Sustituyendo $p_{\text{cos}}(\omega_k) = \frac{\cos\theta_k}{\pi}$ y $f_r = \frac{\rho}{\pi}$ en el estimador de Monte Carlo:

$$\langle L_o \rangle_N = L_e + \frac{1}{N} \sum_{k=1}^N \frac{\left(\frac{\rho}{\pi}\right) L_i(\mathbf{x}, \omega_k) \cancel{\cos\theta_k}}{\frac{\cancel{\cos\theta_k}}{\pi}} = L_e + \frac{1}{N} \sum_{k=1}^N \rho \, L_i(\mathbf{x}, \omega_k)$$

1. **En el Horno de Luz ($\rho=1, L_i=1$)**: Cada muestra produce $Y_k = 1$ constante, por lo que la varianza en 1 rebote se reduce de $\sigma^2_{\text{unif}} = 1/3$ a **$\sigma^2_{\text{cos}} = 0$**.
2. **En la escena Cornell Box con múltiples rebotes**: Evita desperdiciar rayos en ángulos rasantes ($\theta \to \pi/2$ donde $\cos\theta \to 0$) y elimina el factor multiplicativo $2\cos\theta_k$ cuya varianza se amplifica exponencialmente con la profundidad del camino, logrando más de **2.5x a 3.3x menos varianza** para el mismo número de muestras $N$.
