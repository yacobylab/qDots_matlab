function scan = get_FPGA_data(scan,varargin)
global vi
scan.data.FPGA.freqs=double(vi.GetControlValue('freq array')); 
scan.data.FPGA.data=double(vi.GetControlValue('data array'));
end