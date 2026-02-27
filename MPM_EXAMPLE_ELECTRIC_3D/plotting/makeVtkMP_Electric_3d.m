function makeVtkMP_Electric_3d(mpC, phi, J, mpFileName)
%Write VTK file for 3D electric field material point data
%--------------------------------------------------------------------------
% Author: William Coombs (adapted for 3D electric field problem)
% Description:
% Writes ASCII VTK unstructured grid file with:
%   - SCALARS: phi (electric potential)
%   - VECTORS: J   (current density, 3D)
%--------------------------------------------------------------------------
% Input(s):
% mpC        - material point coordinates (nmp,3)
% phi        - electric potential at MPs (nmp,1)
% J          - current density at MPs (nmp,3)
% mpFileName - output file name
%--------------------------------------------------------------------------
[nmp, nD] = size(mpC);
fid = fopen(mpFileName,'wt');
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'3D Electric field MPM output\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i double\n',nmp);
if nD < 3
    mpC3 = [mpC zeros(nmp, 3-nD)];
else
    mpC3 = mpC;
end
fprintf(fid,'%f %f %f \n',mpC3');
fprintf(fid,'\n');
fprintf(fid,'POINT_DATA %i\n',nmp);

% Electric potential (scalar field)
fprintf(fid,'SCALARS phi FLOAT 1\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',phi);
fprintf(fid,'\n');

% Current density (vector field - already 3D)
if nD < 3
    J3 = [J zeros(nmp, 3-nD)];
else
    J3 = J;
end
fprintf(fid,'VECTORS J FLOAT\n');
fprintf(fid,'%f %f %f\n',J3');
fprintf(fid,'\n');

fclose(fid);
end
