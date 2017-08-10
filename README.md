#  Weak Signal Underwater Communications in the Ultra Low Frequency Band in GNU Radio.

This project is related to the paper:

Michel Barbeau, Weak Signal Underwater Communications in the Ultra Low Frequency Band, Proceedings of the GNU Radio Conference, September 2016.


# Copyright 2017 Carleton University.
# Authors: Michel Barbeau
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

Sender: To run a sender, a file in the ".c2" format must be generated first. The program wprsim is used to generate that file , for example:

./wsprsim VE3EMB.c2 "VE3EMB FN25 30"

The source code of wsprsim is available here: 

There are two options. 

![DFLOOD Example](https://github.com/michelbarbeau/gr-splash/blob/master/node_DFLOOD.png)

To run within gnuradio-companion:

Open the flow graph  gr-splash/examples/node_0_DFLOOD.grc

To run outside gnuradio-companion (after generating the flow graph):

cd gr-splash/examples

python top_block.py

Receiver:


Close-loop simulation:





