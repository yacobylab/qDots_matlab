button = questdlg('Ready to quit?', ...
                  'Exit Dialog','Yes','No','No');
         switch button
           case 'Yes',
             disp('Exiting MATLAB');             
           case 'No',
             quit cancel;
end
