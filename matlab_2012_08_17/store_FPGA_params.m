function scan=store_FPGA_params(scan,varargin)
    
global vi;
scan.data.FPGA.Gain=vi.GetControlValue('Gain');
scan.data.FPGA.Threshold=vi.GetControlValue('Threshold control');
scan.data.FPGA.Offset=vi.GetControlValue('Offset');
scan.data.FPGA.freq_indices(1)=vi.GetControlValue('Start 1');
scan.data.FPGA.freq_indices(2)=vi.GetControlValue('End 2');

end