%% Plot compound optic flow vector fields

ml = 3; % mesh limits
[X, Y] = meshgrid(-ml:ml, -ml:ml); % vector field coordinates

% asymmetric optic flow across the eyes
wTrans = 0.5; % weight of translational optic flow
wRot = 1; % weight of yaw optic flow

U_trans = (X./((X.^2 + Y.^2).^0.5))*wTrans; % horizontal vector lengths
V_trans = (Y./((X.^2 + Y.^2).^0.5))*wTrans; % vertical vector lengths

U_trans(4,4) = 0; V_trans(4,4) = 0; % remove NaN

U_yaw = ones(size(X))*wRot;
V_yaw = zeros(size(Y))*wRot;

% plot the vector fields
f = figure('WindowState','maximized'); subplot(1,2,1)
h = quiver(X, Y, U_yaw+U_trans, V_yaw+V_trans);
axis equal; axis off; set(h,'LineWidth',5,'Color','k')

% largely symmetric optic flow across the eyes
wTrans = 5; % weight of translational optic flow
wRot = 1; % weight of yaw optic flow

U_trans = (X./((X.^2 + Y.^2).^0.5))*wTrans; % horizontal vector lengths
V_trans = (Y./((X.^2 + Y.^2).^0.5))*wTrans; % vertical vector lengths

U_trans(4,4) = 0; V_trans(4,4) = 0; % remove NaN

U_yaw = ones(size(X))*wRot;
V_yaw = zeros(size(Y))*wRot;

% plot the vector fields
subplot(1,2,2)
h = quiver(X, Y, U_yaw+U_trans, V_yaw+V_trans);
axis equal; axis off; set(h,'LineWidth',5,'Color','k')