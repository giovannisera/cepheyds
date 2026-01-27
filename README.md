# Hubble Constant Estimation via Cepheid Variables

Questo progetto implementa una pipeline numerica in **Fortran 90** per la stima della **Costante di Hubble ($H_0$)** e dell'età dell'Universo. L'analisi parte dalle curve di luce di 30 stelle Cefeidi osservate in 10 diverse galassie, applicando metodi di calcolo numerico avanzati per la determinazione delle distanze cosmiche.

Il progetto include anche un'analisi cosmologica basata sul modello $\Lambda$-CDM, confrontando i risultati ottenuti con i dati del satellite Planck (2018).

## 🔭 Panoramica del Progetto

Lo scopo è calcolare la velocità di espansione dell'Universo ($H_0$) attraverso la "scala delle distanze", seguendo questi step:
1.  **Calibrazione:** Fit della relazione Periodo-Luminosità (Legge di Leavitt).
2.  **Analisi Curve di Luce:** Determinazione del periodo pulsazionale e magnitudine media tramite interpolazione Spline.
3.  **Distanze:** Calcolo dei moduli di distanza per 10 galassie target.
4.  **Cosmologia:** Stima di $H_0$ e studio dei parametri di densità ($\Omega_m, \Omega_\Lambda$).

**Risultato ottenuto:** $H_0 = 76.54 \pm 1.22 \text{ km s}^{-1} \text{ Mpc}^{-1}$

## 💻 Stack Tecnologico

* **Linguaggio:** Fortran 90
* **Visualizzazione:** Python (Matplotlib)
* **Compilatore:** gfortran (GNU Fortran)

## 🧮 Metodi Numerici Implementati

Il codice è scritto in modo modulare e implementa da zero i seguenti algoritmi:

* **Interpolazione:** Spline Cubica per ricostruire le curve di luce.
* **Fitting:** Minimi quadrati lineari con risoluzione tramite eliminazione di **Gauss-Jordan**.
* **Ottimizzazione:** Minimizzazione del $\chi^2$ ridotto per la determinazione del periodo.
* **Integrazione:** Quadratura di Gauss (Gauss-Legendre a 4 punti) per il calcolo dell'età dell'Universo.
* **Statistica:** Media pesata con propagazione degli errori analitica.

## 📂 Struttura della Repository

* `main.f90`: Programma principale.
* `mod_strumenti.f90`: Modulo strumenti (Sorting, Media Pesata, Spline, Gauss-Jordan).
* `mod_fase1.f90`: Calibrazione relazione P-L.
* `mod_fase2.f90`: Analisi Cefeidi.
* `mod_fase_cosmologia.f90`: Calcolo $H_0$ e integrazione $\Lambda$-CDM.
* `data/`: Cataloghi stellari e curve di luce.

## 🚀 Compilazione ed Esecuzione

Prerequisiti: `gfortran`.

```bash
# Compilazione dei moduli e del main
gfortran -c mod_strumenti.f90 mod_fase1.f90 mod_fase2.f90 mod_fase_cosmologia.f90
gfortran main.f90 mod_strumenti.o mod_fase1.o mod_fase2.o mod_fase_cosmologia.o -o cepheids_analysis

# Esecuzione
./cepheids_analysis
