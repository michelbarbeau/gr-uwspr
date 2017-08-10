#  Weak Signal Underwater Communications in the Ultra Low Frequency Band in GNU Radio.

This project is related to the paper:

Michel Barbeau, Weak Signal Underwater Communications in the Ultra Low Frequency Band, Proceedings of the GNU Radio Conference, September 2017.


# Copyright 2017 Carleton University.
# Author: Michel Barbeau
# Version: August 8, 2017

## Installing 

`git clone https://github.com/michelbarbeau/gr-uwspr`

## Building


```
cd gr-uwspr

mkdir build

cd build 

cmake ../

make

sudo make install

sudo ldconfig

```

## Running

Sender: To run a sender, a file in the ".c2" format must be generated first. The program "wprsim" is used to generate that file, for example:

./wsprsim VE3EMB.c2 "VE3EMB FN25 30"

This program is available in directory "examples".

The source code of "wsprsim" is available here: http://physics.princeton.edu/pulsar/K1JT/wspr.html

The program "wsprsim" generates an audio file in the ".c2" format, named "VE3EMB.c2" in this example. It has in-phase and quadrature signals.

For playing that file, there are two options. 

(1) Within GNU Radio companion, load and run the flowgraph "examples/c2ToAudioSink.grc".

![c2ToAudioSink Example](https://github.com/michelbarbeau/gr-uwspr/blob/master/examples/c2ToAudioSink.png)

(2) Within GNU Radio companion, load and run the flowgraph "examples/c2ToWaveFile".

![c2ToWaveFile Example](https://github.com/michelbarbeau/gr-uwspr/blob/master/examples/c2ToWaveFile.png)

This program converts the ".c2" format into a audio file in the ".wav" format, at 12,000 samples per second. The baseband signal (in the ".c2" file ) is translated to center frequency 1,500 Hz. The output file "test.wav" can be played with any audio player. There is no need to have GNU Radio on the system on which it is played.

Receiver: To run a receiver, load and run the flowgraph "examples/AudioSourceDecode.grc".

![AudioSourceDecode Example](https://github.com/michelbarbeau/gr-uwspr/blob/master/examples/AudioSourceDecode.png)

Closed-loop simulation: A closed-loop simulation is implemented in the flowgraph "examples/WaveFilePlusNoiseDecode.grc". It connects a sender and a noise source to a receiver.

![WaveFilePlusNoiseDecode Example](https://github.com/michelbarbeau/gr-uwspr/blob/master/examples/WaveFilePlusNoiseDecode.png)



