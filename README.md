# SeriesTemporalesIBEX

Código en R usado para realizar el marco práctico del Trabajo de Fin
de Grado **Estudio de la influencia de acciones bancarias y de empresas
nacionales en el IBEX 35 empleando metodologías de series temporales**.

El estudio analiza los log-rendimientos diarios del IBEX 35 y de los cinco
bancos cotizados del selectivo (Banco Santander, BBVA, CaixaBank, Banco
Sabadell y Bankinter) durante el periodo 2015-2025, combinando la metodología
de Box-Jenkins, los modelos ARCH/GARCH de volatilidad condicional, el
análisis cruzado mediante función de transferencia y la cointegración
bilateral y multivariante.

## Estructura del repositorio

El código está organizado en tres scripts independientes que se ejecutan
en orden secuencial.

- `01_univariante.R` : Análisis univariante ARMA + GARCH(1,1) con
  distribución t-Student sobre los seis activos. Para cada activo realiza la
  descarga del precio, el cálculo de log-rendimientos, la limpieza de
  atípicos aditivos, la identificación de la media, el diagnóstico de los
  residuos, la estimación conjunta del modelo, la validación sobre los
  residuos estandarizados y la predicción a corto plazo del precio
  reconstruido con bandas de confianza.

- `02_cruzado.R` : Análisis cruzado en tres bloques, el bilateral por par
  banco vs IBEX 35 (preblanqueo, función de correlación cruzada,
  causalidad de Granger, ajuste Box-Tiao y batería de cointegración
  Engle-Granger, Phillips-Ouliaris y Johansen), el multivariante mediante
  Johansen sobre las seis series, y la función de transferencia agregada con
  los cinco bancos como inputs explicando contemporáneamente al IBEX.

- `03_asimetricos.R` : Bloque opcional con modelización GARCH no lineal
  (sGARCH, GJR-GARCH y EGARCH) restringido a los dos activos donde el Sign
  Bias Test rechaza la simetría del sGARCH(1,1).

## Datos

Las series de precios diarios se descargan en tiempo de ejecución desde
Yahoo Finance vía el paquete `quantmod`. Los tickers utilizados son
`^IBEX`, `SAN.MC`, `BBVA.MC`, `CABK.MC`, `SAB.MC` y `BKT.MC`.


## Reproducibilidad

Los tres scripts fijan `set.seed(42)` al inicio. Los resultados son
estrictamente reproducibles siempre que los datos descargados de Yahoo
Finance coincidan con los del periodo de estudio del TFG (2015-01-01 a
2025-12-31).

## Uso de inteligencia artificial

En cumplimiento de la declaración responsable sobre el uso de inteligencia
artificial, se hace constar que durante el desarrollo de este repositorio se
han empleado herramientas de IA generativa como apoyo en las siguientes tareas:

- **Desarrollo y estructura del código R**: asistencia en la escritura y
  organización de los scripts, en la selección de funciones y paquetes y en la
  estructuración del flujo de análisis.
- **Gráficas y visualizaciones**: asistencia en la generación y el ajuste del
  código encargado de producir las figuras del trabajo.

En todos los casos la herramienta se ha utilizado como apoyo a la
implementación. El diseño del estudio, las decisiones metodológicas, la
elección de los modelos y contrastes estadísticos y la interpretación de los
resultados son responsabilidad del autor, que ha revisado y validado todo el
código y las salidas generadas.
