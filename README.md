# PiperPiTry
This repo contains code samples from my dissertation work in Neuroscience at NYU, completed September 2020. It includes a copy of my dissertation. Within "FINAL_DISSERTATION -- Walter Piper.pdf", Chapter 4 starts on page 77 and the corresponding figures (Figures 4.1-4.22) begin on page 148.

The code samples are either for control of lab equipment (Arduino .ino file or Python .py script) or for data analysis (Python .py script and Jupyter notebooks).

I studied Pavlovian threat conditioning (a form of Classical Conditioning that occurs when auditory tones and mild footshocks are paired). Several generous collaborators connected me to novel biotech: deep-brain fiberoptic recordings (fiber photometry) in rodents expressing a fluorescent norepinephrine receptor. {Special thanks to Yulong Li and his lab at Peking University for sharing the fluorescent receptor technology http://www.yulonglilab.org/publications.html}

This schematic illustrates the fiber photometry set-up.


![Fiber Photometry system](/Figure-fiber-photometry.png)


Extracellular norepinephrine levels in the amygdala were estimated by inputting green fluorescent signal (from the GFP-conjugated norepinephrine receptor) into a moving average model with a control red fluorophore as an exogenous regressor (ARIMAX model of order 0,1,1). All analyses were done with Python. The package statsmodels was used for ARIMAX modeling. Pandas, numpy, and scipy were used for signal processing. Matplotlib was used for visualization.

The plots below show the model trained on 60 seconds of data (turqouise in figure below) and then forecasted up to 30 seconds in the future (orange). The actual future signal (green) was then visualized along with the forecast (orange) (orange-shaded area is 90% confidence interval). The four subplots are consecutive time windows around a single trial of tone-footshock pairing.


![Norepinephrine estimation](/Figure_MovingAvg_ExogRegr_ToneShockPairing.png)


Convolutional neural networks were used to process video data in order to track animal pose, infrared cues, and mechanical stress on the fiberoptic cord near the connection to the fiberoptic implant. This was done with the package DeepLabCut (https://github.com/DeepLabCut/DeepLabCut), which by default uses ResNet-50 with transfer learning for supervised learning of human-labeled anatomical or mechanical features.  Feel free to navigate to DeepLabCut's repo for visualized examples of this deep learning tool in use on various animals.

In the next figure, distance (in pixels) between various points on the rodent's head were plotted on the same time-axis as the fiber photometry signals. 6 distance measurements per frame are plotted below as a time-series scatter plot across the upper section of the figure. Some of the extreme outliers are likely in error, as the (x,y) coordinate estimates were associated with lower Bayesian likelihoods (not shown). The timing of tone delivery (without shock) is indicated by square waves, and the red and green channels of photometry are shown near the top of these square waves, although the plot range on the y-axis is non-optimal for the photometry measures.


![Movement and photometry in full memory test session](/Figure4-17.PNG) 


As an aside, DeepLabCut can also be used to gather (x,y) coordinates from each frame for camera calibration computations and corrections. The figure below shows an example of using the rectilinear geometry of the test cage with monte carlo simulations to estimate the camera calibration parameters.


![Camera correction test](/Figure4-3c.png)


Based on this research, I drew 3 conclusions about the nature of norepinephrine reactivity to psychological stimulation.

Starting from the hypothesis with the strongest confirmatory evidence: Footshocks trigger norepinephrine surges in amygdala (Shown in earlier figures and below). This was previously known based on microdialysis evidence (temporal resolution {15 minute measurements ~ 0.001 Hz}) from the McGaugh Lab in the 1990's, but the research here demonstrated the same with high temporal resolution > 12,000 Hz.

There was also robust evidence for a second hypothesis: Auditory tones, whether already conditioned or novel, tend to elicit robust norepinephrine responses in amygdala. The figure below shows the ARIMAX(0,1,1) output that displays a moderate norepinephrine increase at the start of the 20-second novel tone, a gradual rise during those 20 seconds, and finally a sharp increase during the 1-second footshock that co-terminates with the 20-second tone. The lower plot shows norepinephrine increase during a tone without a shock, then a gradual decline.


![ARIMAX output CSUS and CSminus](/Figure4-16.PNG) 


Finally, there may be a difference between neutral tones and conditionined tones in the maintenance or decay of norepinephrine reactions to successive trials. Norepinephrine reactions to conditioned tones may persist for more trials than neutral tones. This is consistent with theories of "conditioned arousal" or "truncated (internalized) Pavlovian responses".



Thank you for taking the time to read about my past work. If you'd like to get in touch, feel free to add me on LinkedIn (https://www.linkedin.com/in/walter-piper-a320805a/) or message me on ResearchGate (https://www.researchgate.net/profile/Walter_Piper2).

My dissertation lab website: http://www.cns.nyu.edu/ledoux/members.htm
