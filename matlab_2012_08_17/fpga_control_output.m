function []= fpga_control_output(varargin )
%Triggers the fpga to output the control output.
global vi;

if isempty(vi)
    error('vi not valid global variable');
end

vi.SetControlValue('control only',true);
mpause(1e-3);
vi.SetControlValue('Trigger output', true);
mpause(1e-3);
vi.SetControlValue('Trigger output', true);
mpause(1e-3);
vi.SetControlValue('control only', false);

end

