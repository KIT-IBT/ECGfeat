# ECGfeed
A collection of ECG feature extraction algorithms for MATLAB.

## Features
| ID  | Feature                                                    |
|-----|------------------------------------------------------------|
| F1  | first statistical moment of the T wave distribution        |
| F2  | second statistical moment of the T wave distribution       |
| F3  | third statistical moment of the T wave distribution        |
| F4  | fourth statistical moment of the T wave distribution       |
| F5  | RT distance (R peak to T peak)                             |
| F6  | RT mid distance - in the case of biphasic waves            |
| F7  | peakedness of the T wave                                   |
| F8  | T wave amplitude                                           |
| F9  | slope of the ascending part of the T wave                  |
| F10 | slope of the descending part of the T wave                 |
| F11 | ratio of first half T wave energy and whole T wave energy  |
| F12 | ratio of second half T wave energy and whole T wave energy |
| F13 | R peak amplitude                                           |
| F14 | R peak energy                                              |
| F15 | ratio R peak energy and R peak amplitude                   |
| F16 | ST segment change (elevation or depression)                |
| F17 | flag for biphasic T waves (0: monophasic, 1: biphasic)     |
| F18 | R peak area under curve                                    |
## Filtering recommendations
Highpass and lowpass filtering influences the morphology of the ECG. This is why this influence was evaluated for the proposed features and recommendations are given to prevent a distortion by wrong filtering.

|                       |  F1  |  F2  |  F3  |  F4  |  F5  |  F6  |  F7  |  F8  |  F9  |
|-----------------------|:----:|:----:|:----:|:----:|:----:|:----:|:----:|:----:|:----:|
| Lowpass cutoff freq.  |  20  |  40  |  40  |  40  |  20  |  20  |  40  |  40  |  40  |
| Highpass cutoff freq. | 0.05 | 0.10 | 0.10 | 0.10 | 0.05 | 0.05 | 0.10 | 0.10 | 0.10 |

|                       |  F10 |  F11 |  F12 |  F13 |  F14 |  F15 |  F16 |  F17 |  F18 |
|-----------------------|:----:|:----:|:----:|:----:|:----:|:----:|:----:|:----:|:----:|
| Lowpass cutoff freq.  |  40  |  40  |  60  |  70  |  50  |  60  |  40  |   -  |  50  |
| Highpass cutoff freq. | 0.10 | 0.10 | 0.30 | 0.40 | 0.20 | 0.30 | 0.10 |   -  | 0.20 |

## Structure
The structure of the repository is as follows: <br>
**./algorithms** contains the feature extraction algorithms <br>
**./dependencies** contains other projects this one is using <br>
**./examples/study** contains a robustness study of the feature algorithms as well as the calculation of recommendations for filtering without interfering the results from the feature extraction <br>
**./examples/patient_data** contains an example showing how a workflow with a clinical signal from [1][2] could look like <br>

## Next releases
Further algorithms can be found in the development branch and will be added to the main branch in one of the next releases.


## Sources
[1] R. Bousseljot, D. Kreiseler, and A. Schnabel, “Nutzung der EKG-Signaldatenbank CARDIODAT der PTB über das Internet,” Biomedizinische Technik/Biomedical Engineering, pp. 317–318, Jan. 2009.<br>
[2] A. L. Goldberger, L. A. Amaral, L. Glass, J. M. Hausdorff, P. C. Ivanov, R. G. Mark, J. E. Mietus, G. B. Moody, C. K. Peng, and H. E. Stanley, “PhysioBank, PhysioToolkit, and PhysioNet: components of a new research resource for complex physiologic signals.,” Circulation, vol. 101, no. 23, pp. E215–20, Jun. 2000.
