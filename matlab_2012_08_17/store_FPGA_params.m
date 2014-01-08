function scan=store_FPGA_params(scan,varargin)
    
global vi;
scan.data.FPGA.Gain=vi.GetControlValue('Gain');
scan.data.FPGA.Threshold=vi.GetControlValue('Threshold control');
scan.data.FPGA.Offset=vi.GetControlValue('Offset');

end