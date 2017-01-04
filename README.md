{\rtf1\ansi\ansicpg1252\cocoartf1348\cocoasubrtf170
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural

\f0\b\fs24 \cf0 Forward Error Correction in Underwater Acoustic Sensor Networks
\b0 \
\
The project\'92s goal is to analyze the performance of Forward Error Correction (FEC) in Underwater Acoustic Sensor Networks (UASNs). FEC is integrated in GNU Radio environment and tested using Arctic-like simulation. The evaluation is performed according to the message size after encoding and Packet Error Rate (PER) vs. Signal to Noise Ratio (SNR).\
\

\b Brief instruction:\
\

\b0 The project requires GNU Radio and Matlab installed.\
\
1. Start GNU Radio companion. Open one of the FSK Sender .grc files. The project contains three versions of sender and receiver: FSK Sender/Receiver no FEC, FSK Sender/Receiver with FEC Tagged CC, FSK Sender/Receiver with Asynchronous POLAR. \
\
2. In FSK Sender specify the source file with data.\
\
3. Start the file. The file generates FSK_Output.wav file.\
\
4. Open Matlab. Add FSK_Output.wav to 5km_80pc_ice folder. Add a path to 5km_80pc_ice and at folders. Run pskProcessing command with the number of the environment file, the duration of wav file in seconds, and \'91GFSK\'92 as parameters.\
\
5. Matlab generates 17 wav files for each environment file. Open GNU Radio and FSK Receiver file in gnu radio. In a wav file source specify a path to one of the wav files and run.\
\
More details are provided in my honours paper.}