#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Fsk Sender
# Generated: Wed Jul  6 01:51:58 2016
##################################################

from gnuradio import blocks
from gnuradio import digital
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import math


class FSK_Sender(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Fsk Sender")

        ##################################################
        # Variables
        ##################################################
        self.transistion = transistion = 100
        self.sps = sps = 8
        self.sideband_rx = sideband_rx = 1000
        self.sideband = sideband = 1000
        self.samp_rate = samp_rate = 48000
        self.interpolation = interpolation = 18000
        self.decimation = decimation = 1
        self.carrier = carrier = 10000

        ##################################################
        # Blocks
        ##################################################
        self.rational_resampler_xxx_0_1 = filter.rational_resampler_ccc(
                interpolation=interpolation,
                decimation=decimation,
                taps=None,
                fractional_bw=None,
        )
        self.freq_xlating_fir_filter_xxx_0 = filter.freq_xlating_fir_filter_ccc(1, (firdes.band_pass (0.5,samp_rate,10000-sideband,10000+sideband,transistion)), -carrier, samp_rate)
        self.digital_gfsk_mod_0 = digital.gfsk_mod(
        	samples_per_symbol=sps,
        	sensitivity=1.0,
        	bt=0.35,
        	verbose=False,
        	log=False,
        )
        self.blocks_wavfile_sink_1 = blocks.wavfile_sink("FSK_output.wav", 1, 48000, 16)
        self.blocks_unpack_k_bits_bb_0_0 = blocks.unpack_k_bits_bb(8)
        self.blocks_throttle_0 = blocks.throttle(gr.sizeof_gr_complex*1, samp_rate,True)
        self.blocks_multiply_const_vxx_0_0 = blocks.multiply_const_vff((0.1, ))
        self.blocks_file_source_0 = blocks.file_source(gr.sizeof_char*1, "TestData2", False)
        self.blocks_file_sink_0_0 = blocks.file_sink(gr.sizeof_char*1, "inputBinary", False)
        self.blocks_file_sink_0_0.set_unbuffered(False)
        self.blocks_complex_to_real_0_0 = blocks.complex_to_real(1)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_complex_to_real_0_0, 0), (self.blocks_multiply_const_vxx_0_0, 0))    
        self.connect((self.blocks_file_source_0, 0), (self.blocks_unpack_k_bits_bb_0_0, 0))    
        self.connect((self.blocks_file_source_0, 0), (self.digital_gfsk_mod_0, 0))    
        self.connect((self.blocks_multiply_const_vxx_0_0, 0), (self.blocks_wavfile_sink_1, 0))    
        self.connect((self.blocks_throttle_0, 0), (self.rational_resampler_xxx_0_1, 0))    
        self.connect((self.blocks_unpack_k_bits_bb_0_0, 0), (self.blocks_file_sink_0_0, 0))    
        self.connect((self.digital_gfsk_mod_0, 0), (self.blocks_throttle_0, 0))    
        self.connect((self.freq_xlating_fir_filter_xxx_0, 0), (self.blocks_complex_to_real_0_0, 0))    
        self.connect((self.rational_resampler_xxx_0_1, 0), (self.freq_xlating_fir_filter_xxx_0, 0))    

    def get_transistion(self):
        return self.transistion

    def set_transistion(self, transistion):
        self.transistion = transistion
        self.freq_xlating_fir_filter_xxx_0.set_taps((firdes.band_pass (0.5,self.samp_rate,10000-self.sideband,10000+self.sideband,self.transistion)))

    def get_sps(self):
        return self.sps

    def set_sps(self, sps):
        self.sps = sps

    def get_sideband_rx(self):
        return self.sideband_rx

    def set_sideband_rx(self, sideband_rx):
        self.sideband_rx = sideband_rx

    def get_sideband(self):
        return self.sideband

    def set_sideband(self, sideband):
        self.sideband = sideband
        self.freq_xlating_fir_filter_xxx_0.set_taps((firdes.band_pass (0.5,self.samp_rate,10000-self.sideband,10000+self.sideband,self.transistion)))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.blocks_throttle_0.set_sample_rate(self.samp_rate)
        self.freq_xlating_fir_filter_xxx_0.set_taps((firdes.band_pass (0.5,self.samp_rate,10000-self.sideband,10000+self.sideband,self.transistion)))

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
        self.freq_xlating_fir_filter_xxx_0.set_center_freq(-self.carrier)


def main(top_block_cls=FSK_Sender, options=None):

    tb = top_block_cls()
    tb.start()
    tb.wait()


if __name__ == '__main__':
    main()
