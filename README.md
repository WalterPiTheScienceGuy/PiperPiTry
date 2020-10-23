# PiperPiTry
This repo contains code samples from my dissertation work in neural science. It also contains a copy of my dissertation (Figures for fiber photometry begin on page 136/148).

I studied Pavlovian threat conditioning (a form of Classical Conditioning). the project as a whole involved fiberoptic recording (fiber photometry) from non-human subjects expressing a fluorescent norepinephrine receptor. The code samples are either for control of lab equipment (Arduino or Python) or for data analysis (Python script and Jupyter notebook).

![Fiber Photometry system](/Figure-fiber-photometry.png)

Extracellular norepinephrine levels in the amygdala were estimated by inputting green fluorescent signal (from the GFP-conjugated norepinephrine receptor) into a moving average model with a control red fluorophore as an exogenous regressor (ARIMAX model of order 0,1,1). The Python package statsmodels was used for the ARIMAX modeling.

The plots below show the model trained on 60 seconds of data (turqouise in figure below) and then forecasted up to 30 seconds in the future (orange). The actual future signal (green) was then visualized along with the forecast (orange) (orange-shaded area is 90% confidence interval). The four subplots are consecutive time windows around a single experimental tone-footshock pairing.

![Norepinephrine estimation](/Figure_MovingAvg_ExogRegr_ToneShockPairing.png)

Convolutional neural networks were used to process video data in order to track animal pose, infrared cues, and mechanical stress on the fiberoptic cord near the connection to the fiberoptic implant.
