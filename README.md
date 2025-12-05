=====================================================================
                           K D _ F A S T
=====================================================================

---------------------------------------------------------------------
PROGRAM PARAMETER DESCRIPTION
---------------------------------------------------------------------

The program supports several optional parameters to customize the
hydrophobicity profile calculation (Kyte–Doolittle or any
sliding-window–based method). 
Time for 750.000 protein seq (i7 threads)

real	0m4.460s
user	0m34.574s
sys	0m0.589s

multi-thread high-performance processing model described earlier.

Available parameters:

./kd_fast_openmp [input] [cutoff] [win_min=19] [win_max=21] [threads]
Example:
./kd_fast_openmp all_idx 1.6 19 21 8


    [cutoff]
    [win_min=19]
    [win_max=21]
    [threads]

A detailed explanation of each parameter follows.

OUTPUT summary:
 - Column 1: sequence ID
 - Column 2: hydropathy values for each window position
 - Column 3: left position
 - Column 4: right position
 - Column 5: windows slide length
 
---------------------------------------------------------------------
[cutoff]
---------------------------------------------------------------------
Defines a numerical threshold used to filter or highlight hydrophobic
regions in the computed profile. Typical uses include:

 - Detecting long hydrophobic stretches
 - Identifying potential transmembrane helices
 - Highlighting peaks above biologically relevant levels
 - Emitting only windows whose mean score exceeds the threshold

A common Kyte–Doolittle threshold:
    cutoff ≈ 1.6
(Used in many membrane-protein predictions)

When `cutoff` is provided, the program may:
 - Mark positions above the threshold
 - Output only windows that exceed the value
 - Produce a binary profile (0/1)
 - Suppress non-significant windows

If no cutoff is given, the full unfiltered profile is printed.

---------------------------------------------------------------------

[win_min = 19]
[win_max = 21]

---------------------------------------------------------------------
These parameters define the **range of sliding-window sizes** to be
computed for each sequence.

Although the original KD method used windows of size 9, 11, or 19,
modern analyses of membrane proteins commonly use:

 - Long windows → 19–21 amino acids
 - Smoother profiles with reduced noise
 - Enhanced detection of broad hydrophobic peaks
 - Better discrimination of true transmembrane segments

The program accepts:
 - `win_min` → minimum window size
 - `win_max` → maximum window size

If a range is provided (e.g., 19–21):
    The program computes all window lengths in that range:
        19 aa
        20 aa
        21 aa

This comparative multi-window approach helps identify stable
hydrophobic segments independent of the exact window length.

If a fixed size is desired:
    win_min = win_max = <value>

---------------------------------------------------------------------
[threads]
---------------------------------------------------------------------
Specifies the number of CPU threads to use during processing.  
The program is designed for high-throughput workloads and scales well
on multi-core processors.

 - threads = 1  
       → sequential execution
 - threads > 1  
       → parallel processing of independent sequences

Advantages:
 - Near-linear speedup up to ~8 threads (typical for older i7 CPUs)
 - Strong performance for very large datasets (tens of millions of seqs)
 - Fully utilizes modern and legacy multi-core systems

Default behavior (implementation-dependent):
 - If not specified: use (number_of_cores - 1)
 - `--threads auto` may automatically detect the optimal value

Typical usage:
    --threads 8
    --threads 4
    --threads auto

---------------------------------------------------------------------
QUICK SUMMARY
---------------------------------------------------------------------

cutoff    : threshold for hydrophobic windows (KD typical >1.6)
win_min   : minimum sliding-window size (19 recommended for TM regions)
win_max   : maximum sliding-window size (21 recommended for TM regions)
threads   : number of CPU threads used for parallel processing


---------------------------------------------------------------------
1. DESCRIPCIÓN GENERAL
---------------------------------------------------------------------
KD_FAST es una herramienta ultrarrápida diseñada para calcular
perfiles de hidrofobicidad (Kyte-Doolittle basados en
ventanas deslizantes) sobre millones de secuencias de proteínas.

El software prioriza:
 - Velocidad extrema
 - Bajo overhead
 - Uso de memoria optimizado

Soporta secuencias delimitadas por espacio:
   seq_id  protein_sequence

Ejemplo:
   P12345  MAVLVLGLAALAT...

---------------------------------------------------------------------
2. ESTIMACIÓN DE RENDIMIENTO (calculada por GPT)
---------------------------------------------------------------------
Hardware estimado:
 - Intel Core i7 (8 threads) modelo antiguo (2015–2017)
 - 16 GB RAM sugeridos para 100 millones de secuencias

Carga de trabajo:
 - 100 millones de secuencias
 - Longitud promedio = 1000 aa
 - Ventana = 19 (Kyte-Doolittle típica)

Cómputo estimado:
 - ~100 millones × (1000 - 19) ≈ 9.81 × 10^10 posiciones analizadas
 - i7 antiguo: 280–350 millones de operaciones/seg por thread práctico
 - Escalamiento a 8 threads: 2.0–2.5 × 10^9 ops/s

Tiempo estimado:
 - ~98 × 10^9 operaciones / 2.2 × 10^9 ops/s ≈ 44–52 segundos

Consumo esperado:
 - RAM usada principalmente por buffers de entrada: ~1–2 GB
 - CPU al 100% en 8 threads

(*Todos estos cálculos fueron generados automáticamente por ChatGPT a
solicitud del usuario y dependen de implementación y compilación.*)

Para 717.615 secuencia de proteinas de tamaño variable:

Tiempo a 8 hilos de CPU:
real	0m4.460s
user	0m34.574s
sys	0m0.589s

---------------------------------------------------------------------
3. FUNCIONAMIENTO DE LA VENTANA DESLIZANTE (SLIDING WINDOW)
---------------------------------------------------------------------

Una "ventana" es un segmento de tamaño fijo (ej. 19 aa) que se mueve
sobre la secuencia calculando un valor por cada posición.
Esto produce un perfil continuo que permite detectar regiones
hidrofóbicas, dominios transmembrana, señales, picos, valles, etc.

---------------------------------------------------------------------
4. RESUMEN EN ASCII + GRÁFICO EXPLICATIVO (TOTALMENTE EN ASCII)
---------------------------------------------------------------------

RESUMEN:
---------
La ventana deslizante recorre una secuencia tomando subsegmentos
consecutivos. En cada paso calcula un valor (como hidrofobicidad).
Esto suaviza el ruido y revela patrones locales o regiones funcionales.

GRÁFICO ASCII DE VENTANA DESLIZANTE:
------------------------------------

SECUENCIA EJEMPLO:
A  R  L  V  I  G  D  T  Y  W  S  E  M  K  F
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |

VENTANA = 5

Posición 1:
┌─────────────────┐
│ A  R  L  V  I    │  → valor = [■■■■■■■■░░]
└─────────────────┘
        ↓

Posición 2:
    ┌─────────────────┐
    │ R  L  V  I  G    │ → valor = [■■■■■■■■■■]
    └─────────────────┘
           ↓

Posición 3:
       ┌─────────────────┐
       │ L  V  I  G  D    │ → valor = [■■■■■░░░░]
       └─────────────────┘
              ↓

Posición 4:
          ┌─────────────────┐
          │ V  I  G  D  T    │ → valor = [■■■■░░░░░]
          └─────────────────┘
                 ↓

---------------------------------------------------------------------
8. LICENCIA
---------------------------------------------------------------------
Software abierto para uso académico e industrial.
El usuario es libre de modificar, mejorar o integrar KD_FAST en
cualquier pipeline de bioinformática.
Posición 5:
             ┌─────────────────┐
             │ I  G  D  T  Y    │ → valor = [■■■■■■■░░]
             └─────────────────┘

Representación global:
A R L V I G D T Y W S E M K F
[#####-----]
  [######----]
    [#####-----]
      [###-------]
        [####------]
          ...

Esto ilustra cómo KD_FAST evalúa cada ventana para generar
perfiles grandes, suaves y altamente informativos.

---------------------------------------------------------------------
5. COMPILACIÓN E INSTALACIÓN
---------------------------------------------------------------------

Usando make:
   $ make

Usando gcc directamente:
   $ gcc -O3 -march=native -mtune=native -o kd_fast kd_fast.c

Instalar en /usr/local/bin:
   $ sudo cp kd_fast /usr/local/bin/

---------------------------------------------------------------------
6. ENTRADA Y SALIDA
---------------------------------------------------------------------

ENTRADA (STDIN):
   seq_id sequence

SALIDA (STDOUT):
   seq_id valor1 valor2 valor3 ...

MANEJO DE ERRORES:
 - Si la línea no contiene exactamente dos campos, se imprime:
   ERROR: formato inválido en la línea X

---------------------------------------------------------------------
7. NOTAS DE IMPLEMENTACIÓN
---------------------------------------------------------------------
 - El código evita malloc dentro del loop crítico.
 - Se usa tabla estática de hidrofobicidad.
 - Ciclo ajustado para máxima velocidad.
 - Se recomienda compilar con GCC ≥ 9 para mejor vectorización.
 - Se puede añadir OpenMP para paralelismo masivo.

=====================================================================
                       FIN DEL README
=====================================================================
