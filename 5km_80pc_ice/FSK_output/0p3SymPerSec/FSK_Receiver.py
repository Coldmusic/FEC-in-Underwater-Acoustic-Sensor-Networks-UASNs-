#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Fsk Receiver
# Generated: Wed Jul  6 01:54:12 2016
##################################################
from sys import argv
from gnuradio import blocks
from gnuradio import digital
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import math


class FSK_Receiver(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Fsk Receiver")

        ##################################################
        # Variables
        ##################################################
        self.transistion = transistion = 100
        self.sps = sps = 8
        self.sideband_rx = sideband_rx = 500
        self.sideband = sideband = 500
        self.samp_rate = samp_rate = 48000
        self.interpolation = interpolation = 1
        self.decimation = decimation = 1
        self.carrier = carrier = 10000

        ##################################################
        # Blocks
        ##################################################
        self.rational_resampler_xxx_0_0_0 = filter.rational_resampler_ccc(
                interpolation=decimation,
                decimation=100,
                taps=None,
                fractional_bw=None,
        )
        self.rational_resampler_xxx_0_0 = filter.rational_resampler_ccc(
                interpolation=decimation,
                decimation=180,
                taps=None,
                fractional_bw=None,
        )
        self.freq_xlating_fir_filter_xxx_0_0_0 = filter.freq_xlating_fir_filter_ccc(1, (filter.firdes.low_pass(1, samp_rate*10, sideband_rx,1000)), carrier, samp_rate)
        self.digital_gfsk_demod_0 = digital.gfsk_demod(
        	samples_per_symbol=sps,
        	sensitivity=1.0,
        	gain_mu=0.175,
        	mu=0.5,
        	omega_relative_limit=0.005,
        	freq_error=0.0,
        	verbose=False,
        	log=False,
        )
        
        
        script, inputwav, outputBinary= argv


        self.blocks_wavfile_source_0 = blocks.wavfile_source(inputwav, False)
        self.blocks_throttle_1_0_0 = blocks.throttle(gr.sizeof_gr_complex*1, samp_rate,True)
        self.blocks_float_to_complex_0 = blocks.float_to_complex(1)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_char*1, outputBinary, False)
        self.blocks_file_sink_0.set_unbuffered(False)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_float_to_complex_0, 0), (self.blocks_throttle_1_0_0, 0))    
        self.connect((self.blocks_throttle_1_0_0, 0), (self.freq_xlating_fir_filter_xxx_0_0_0, 0))    
        self.connect((self.blocks_wavfile_source_0, 0), (self.blocks_float_to_complex_0, 0))    
        self.connect((self.digital_gfsk_demod_0, 0), (self.blocks_file_sink_0, 0))    
        self.connect((self.freq_xlating_fir_filter_xxx_0_0_0, 0), (self.rational_resampler_xxx_0_0, 0))    
        self.connect((self.rational_resampler_xxx_0_0, 0), (self.rational_resampler_xxx_0_0_0, 0))    
        self.connect((self.rational_resampler_xxx_0_0_0, 0), (self.digital_gfsk_demod_0, 0))    

    def get_transistion(self):
        return self.transistion

    def set_transistion(self, transistion):
        self.transistion = transistion

    def get_sps(self):
        return self.sps

    def set_sps(self, sps):
        self.sps = sps

    def get_sideband_rx(self):
        return self.sideband_rx

    def set_sideband_rx(self, sideband_rx):
        self.sideband_rx = sideband_rx
        self.freq_xlating_fir_filter_xxx_0_0_0.set_taps((filter.firdes.low_pass(1, self.samp_rate*10, self.sideband_rx,1000)))

    def get_sideband(self):
        return self.sideband

    def set_sideband(self, sideband):
        self.sideband = sideband

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.blocks_throttle_1_0_0.set_sample_rate(self.samp_rate)
        self.freq_xlating_fir_filter_xxx_0_0_0.set_taps((filter.firdes.low_pass(1, self.samp_rate*10, self.sideband_rx,1000)))

    def get_interpolation(self):
        return self.interpolation

    def set_interpolation(self, interpolation):
        self.interpolation = interpolation

    def get_decimation(self):
        return self.decimation

    def set_decimation(self, decimation):
        self.decimation = decimation

    def get_carrier(self):
        return self.carrier

    def set_carrier(self, carrier):
        self.carrier = carrier
        self.freq_xlating_fir_filter_xxx_0_0_0.set_center_freq(self.carrier)


def main(top_block_cls=FSK_Receiver, options=None):

    tb = top_block_cls()
    tb.start()
    tb.wait()


if __name__ == '__main__':
    main()
