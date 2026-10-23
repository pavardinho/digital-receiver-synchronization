# Digital Receiver Synchronization and Image Reconstruction

An academic MATLAB project implementing a receiver chain for a noisy, oversampled digital signal and reconstructing the transmitted grayscale image.

## Receiver pipeline

1. Apply a square-root raised-cosine matched filter and compensate its group delay.
2. Detect the useful signal interval and select the symbol timing phase from the received energy.
3. Track and correct carrier phase/frequency error with a second-order phase-locked loop (PLL).
4. Demodulate the symbols, decode the image bits, and display the reconstructed image.

`code_recepteur.m` runs the receiver. `tests.m` contains exploratory plots used to inspect the received waveform, spectrum, matched filter, timing, and constellation.

## Input data

The course-provided `signal_recu.mat` file is intentionally not included in this portfolio copy. To run either script, place an authorized copy beside the MATLAB files. It must contain a variable named `signal_recu` holding the received complex signal. The scripts display a clear error if the file or variable is missing.

## Requirements

MATLAB with Communications Toolbox (`rcosdesign`, `bi2de`) is required for the receiver. The analysis script also uses `pwelch` from Signal Processing Toolbox.

Run `code_recepteur.m` in MATLAB to reconstruct and display the image. Run `tests.m` to generate the signal-analysis figures.

## Scope

This repository contains the receiver implementation and analysis scripts only; it excludes course handouts, presentation slides, and received-signal data.
