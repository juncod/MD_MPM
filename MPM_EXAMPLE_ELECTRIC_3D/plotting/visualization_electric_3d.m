function visualization_electric_3d(mpData, mesh)
%Visualization of 3D electric field MPM results
%--------------------------------------------------------------------------
% Plots electric potential distribution and current density vectors
% at material points using 3D scatter plots.
%--------------------------------------------------------------------------
nmp = length(mpData);
pos = zeros(nmp,3);
phi_vals = zeros(nmp,1);
Jx = zeros(nmp,1);
Jy = zeros(nmp,1);
Jz = zeros(nmp,1);

for i=1:nmp
    pos(i,:) = mpData(i).mpC;
    phi_vals(i) = mpData(i).phi;
    Jx(i) = mpData(i).J(1);
    Jy(i) = mpData(i).J(2);
    Jz(i) = mpData(i).J(3);
end

figure('Position',[100 100 1400 500]);

% Left panel: electric potential (scalar field) - 3D scatter
subplot(1,2,1);
scatter3(pos(:,1), pos(:,2), pos(:,3), 8, phi_vals, 'filled');
grid on; axis equal; axis tight;
colorbar; title('\phi  (Electric Potential)');
xlabel('x'); ylabel('y'); zlabel('z');
view(3);

% Right panel: current density magnitude with vector arrows
subplot(1,2,2);
J_mag = sqrt(Jx.^2 + Jy.^2 + Jz.^2);
scatter3(pos(:,1), pos(:,2), pos(:,3), 8, J_mag, 'filled');
grid on; axis equal; axis tight;
colorbar; title('|J|  (Current Density Magnitude)');
xlabel('x'); ylabel('y'); zlabel('z');
hold on;
step = max(1, floor(nmp/200));
quiver3(pos(1:step:end,1), pos(1:step:end,2), pos(1:step:end,3), ...
        Jx(1:step:end), Jy(1:step:end), Jz(1:step:end), 0.5, 'k');
view(3);
end
