# MATLAB code to replicate main results from "Decoding Motor Plans Using a Closed-Loop Ultrasonic Brain-Machine Interface"
==============================


## Project Organization
------------

    ├── README.md           <- The top-level README for developers using this project
    ├── setup.m             <- Used to setup the repo and data directory hierarchy
    ├── behavior            <- functions for parsing behavioral files and data
    ├── data loading        <- functions for loading data
    ├── decoders            <- functions related to decoding, including different decoder algorithms
    │   └── classwise PCA   <- specialized functions for classwise PCA
    ├── images              <- Referenced by `README.md` to explain dialog boxes and user choices.
    ├── image alignment     <- functions for aligning different datasets and pretraining models
    ├── plotting            <- functions for visualization of results and output
    ├── statistics          <- functions for calculating different types of statistics
    ├── utilities           <- functions with miscellanous purposes
    └── third-party         <- third party repositories, functions, etc.

------------
## Associated dataset
Available from CaltechData: [Click here to download](https://doi.org/10.22002/pa710-cdn95). 

Packaged as .zip file in data hierarchy used by MATLAB functions. Once downloaded, extract contents from the .zip file.
   
    ├── project_record.json               Metadata about each experimental session.
    ├── DescriptionOfFiles.pdf            Information about contents of .mat files
    ├── aligned_doppler_data              Empty folder where aligned datasets will be saved
    ├── doppler data
    |   └──rt_fUS_data_S{*A}_R{*B}.mat    Doppler data files where {*A} is the session number and {*B} is the run number.
    ├── output                            Empty folder where figures and output will be saved.
    ├── simulated                         Empty folder where simulated BMI sessions will be saved.
    ├── sulcus                            Empty folder where sulcal maps will be saved.

------------
## Other information
Tested on MATLAB R2021a on Windows 10, MATLAB 2021b, and MATLAB 2022b on Mac Ventura 13.4.1
Please send feedback and suggestions to: [wsgriggs@caltech.edu](mailto:wsgriggs@caltech.edu)

Zenodo archive - 
[![DOI](https://zenodo.org/badge/701341316.svg)](https://zenodo.org/badge/latestdoi/701341316)

==============================

### In publications, please reference:
Griggs, W.S., Norman, S.L., Deffieux, T., Segura, F., Osmanski, B.-F., Chau, G., Christopoulos, V., Liu, C., Tanter, M., Shapiro, M.G., and Andersen, R.A. Decoding Motor Plans Using a Closed-Loop Ultrasonic Brain-Machine Interface. Nature Neuroscience. November 30, 2023. https://www.nature.com/articles/s41593-023-01500-7

------------
## Quick Start
1. Clone this repo to a known location on your computer.

2. Download paired dataset from [CaltechData](https://doi.org/10.22002/pa710-cdn95) and unzip in a location that is convenient.

4. Run `setup.m`. This will add appropriate folders to MATLAB search path and specify path to the downloaded data directory.

5. Run code blocks from  `main_analysis.m`.

------------
## Demo of `simulate_real_time_fUS_BMI`
1. Select fUS-BMI run mode. Here, I chose the verbose mode.
* `Run as fast as possible with minimal user display` - Fast mode; Will display occasional updates to the command line and MATLAB figure about performance
* `Run slower with full display of data streaming and trial performance` - Verbose mode; Will display simulated data streaming, task, and trial-by-trial performance of fUS-BMI
![Dialog box to select fUS-BMI run mode](images/1_fUS_BMI_mode_selection.png)

2. Select training mode. Here, I chose the `pretrain+retrain` option.
* `pretrain+retrain` - Use both a pretrained model and retrain after each trial
* `pretrain only` - Use only a pretrained model and do not retrain model.
* `retrain only` - Use only data from current session and retrain after each trial after sufficient trials for initial training set.
* For more description of pretraining vs retraining, see paper.
  
![Dialog box to select training mode](images/2_select_training_mode.png)

1. Select dataset to test (and retrain) fUS-BMI on.
* If `retrain` option was selected in step #2. Then data from this dataset will also be included in training dataset.
* Does not matter whether the underlying dataset was collected with `pretrain` or `retrain` on. All of the data is compatible with all of the modes.
 
![Dialog box to select test dataset](images/3_dataset_selection.png)

4. Select if you want to save the data from this run
![Dialog box to select if to save data](images/4_save_data.png)

5. If you want to save the data, select where to save.
![Dialog box to select where to save data](images/5_where_to_save_data.png)

6. If you specified to use `pretrain` option, then select aligned dataset to use
* If no aligned datasets exist or you want to create a new one, then select cancel.


![Dialog box to specify aligned dataset](images/6_specify_aligned_dataset.png)

7. Specify whether you want to create an aligned dataset to be used for pretraining or if you just want to quit.
![Dialog box to create new alignment](images/7_create_new_alignment.png)

8. Specify which dataset you want to use for pretraining.
> [!WARNING]
> If you choose the same dataset you chose earlier, this will create a data leak.

![Dialog box to select alignment dataset](images/8_select_alignment_dataset.png)

10. Align neurovascular maps using the MATLAB app.
* `Automated Intensity-Based Image Alignment` button - This should get the two neurovascular images mostly aligned.
* If the alignment can be improved, use the arrow keys to translate the image and `r`/`t` to rotate the image.
* Once you are happy with the alignment, click `save transform to workspace` button.
* If you want to reset to initial alignment, click `reset to initial alignment` button 

### Pre-alignment
![MATLAB App to align datasets](images/9_align_datasets.png)

### Post-alignment
![Post alignment](images/9b_aligned_datasets.png)

10. Select where to save the aligned dataset.
![Dialog box to select where to save aligned data](images/10_where_to_save_aligned_dataset.png)

12. The fUS-BMI will now run.
![fUS-BMI now running](images/11_example_display.png)

------------
## Demo of `simulate_real_time_fUS_BMI`
1. Select dataset
![Select dataset](images/SL_0_select_dataset.png)
   
1.5. (if sulcus has been previously traced and saved) Can specify if you want to create a new sulcus map
![Create new sulcus map?](images/SL_05_define_new_sulcus.png)

2. Trace the brain surface and sulci.
* Trace the brain surface and sulci by sequentially clicking on the image with the crosshairs.
* Once finished, press Enter.
* There should be one line between the left and right side of the image and the line should not cross itself.
  
### Pre-tracing
![Trace new sulcus](images/SL_1_trace_sulcus.png)

### Part way through tracing
![Partway through drawing sulcus](images/SL_1b_trace_sulcus.png)

### Post-tracing
![Finished drawing sulcus](images/SL_1c_trace_sulcus.png)

3. Save the sulcus map to file
![Save sulcus map](images/SL_2_save_sulcus.png)

4. Save the searchlight results to file
![Save searchlight results](images/SL_3_save_results.png)

5. Display of searchlight analysis results
![Display searchlight results](images/SL_4_searchlight_results.png)


