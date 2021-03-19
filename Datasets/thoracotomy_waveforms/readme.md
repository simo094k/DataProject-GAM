# Waveforms recorded before and after thoracotomy

## Windows

closed_chest* -> *Fluid* -> closed_chest(post fluid) -> *Cut open thorax (thoracotomy)* -> open_chest -> *Fluid* -> open_chest(post fluid)*

*not in this dataset.

## Data structure (individual pt waveforms)

```
List of 3
 $ id   : Patient ID
 $ closed_chest: Data recorded with closed chest
  ..$ ECG                : Electrocardiogram (500 Hz)
  ..$ ABP                : Arterial blood pressure (125 Hz)
  ..$ CVP                : Central venous pressure (125 Hz)
  ..$ QRSmarker_ecg_index: Occurrence of QRS complex in ECG vector (index)  
  ..$ QRSmarker_ms       : Occurrence of QRS complex in milliseconds from start
 $ open_chest  : Data recorded with closed chest
  ..$ ECG                : ... Same as above
  ..$ ABP                : 
  ..$ CVP                : 
  ..$ QRSmarker_ecg_index: 
  ..$ QRSmarker_ms       : 
```

## comb_ventilation_data

Data frame of ventilator settings for all patients. This is "validated" for pre_fluid closed chest. It is probably the same post fluid.