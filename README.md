#  Low Frequency and Weak Signal Underwater Communications in GNU Radio.

See the Wiki section of this repository for a list of related papers and MATLAB Live Scripts explaining some of the details of the software design.

# Copyright 2019 Michel Barbeau, Carleton University.
# Author: Michel Barbeau
# Version: February 22, 2019

# Latest build on GNURadio 3.8
# updated by Jay Patel
# version: June 26, 2020.

## Installing

`git clone https://github.com/patel999jay/gr-uwspr.git`

Note : This is just fork of the original repo : https://github.com/michelbarbeau/gr-uwspr . Please refer the same for more details.

## Building


```
cd gr-uwspr

mkdir build

cd build

cmake -DENABLE_DOXYGEN=OFF -DGNURADIO_ALL_INCLUDE_DIRS=/usr/include/gnuradio/swig ../

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

Note: Paths to local files in "Wav File Source" blocks must be updated according to your local installation.

Watch the following demo on Youtube.

[![Demo](https://i1.ytimg.com/vi/98o4X0QdZ78/hqdefault.jpg)](https://youtu.be/98o4X0QdZ78)
