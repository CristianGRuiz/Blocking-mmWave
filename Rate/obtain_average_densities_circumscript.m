function [E_L, lambda_b, NumberBuildings, A_T, F_L_emp, l_emp_cdf] = obtain_average_densities_circumscript()

load('ChicagoDowntownBuildingStats.mat')

x_coord = [];
y_coord = [];

% We loop over all the buildings of the database
for n_building = 1:NumberBuildings
    % We obtain the number of edges that form the building n_building
    N_edges = sum(BuildingXYCoord(n_building,:,1)~=0)-1;
    % Then, we loop over all the edges that form a given building
    for n_edge = 1:N_edges
        % We obtain the 2 values of x_coord that define an edge
        x_coord = [x_coord; BuildingXYCoord(n_building,n_edge,1) BuildingXYCoord(n_building,n_edge+1,1)];
        % We obtain the 2 values of y_coord that define an edge
        y_coord = [y_coord; BuildingXYCoord(n_building,n_edge,2) BuildingXYCoord(n_building,n_edge+1,2)];
    end
end

% We obtain the angles of each segment with the current coordinates system
theta = mod(angle(exp(1i*(angle(x_coord(:,2)+1i*y_coord(:,2)-(x_coord(:,1)+1i*y_coord(:,1)))))),pi);

% We obtain the angles of each segment with the current coordinates system
theta(theta>deg2rad(10),:) = [];

% We obtain the reference angle to which all segments are turned
theta_ref = mean(theta);

% The modulus from the origin to the first point of the segment
mod_1 = sqrt(x_coord(:,1).^2+y_coord(:,1).^2);
mod_2 = sqrt(x_coord(:,2).^2+y_coord(:,2).^2);

% Now, the angles of all the points with respect of the origin are
% computed
alpha_1 = angle(exp(1i*(angle(x_coord(:,1)+1i*y_coord(:,1)))));
alpha_2 = angle(exp(1i*(angle(x_coord(:,2)+1i*y_coord(:,2)))));

% Then, the scenario is rotated
x_coord = [mod_1.*cos(alpha_1-theta_ref), mod_2.*cos(alpha_2-theta_ref)];
y_coord = [mod_1.*sin(alpha_1-theta_ref), mod_2.*sin(alpha_2-theta_ref)];

% figure(1)
% hold on;
% [N_segments,~] = size(x_coord); 
% for n_segment = 1:N_segments
%     plot(x_coord(n_segment,:),y_coord(n_segment,:),'b','LineWidth',2)
% end

% The modulus from the origin to the first point of the segment
norma = sqrt(BuildingXYCoord(:,:,1).^2+BuildingXYCoord(:,:,2).^2);

% Now, the angles of all the points with respect of the origin are
% computed
alpha = angle(exp(1i*(angle(BuildingXYCoord(:,:,1)+1i*BuildingXYCoord(:,:,2)))));

% Then, the scenario is rotated
BuildingXYCoord(:,:,1) = norma.*cos(alpha-theta_ref);
BuildingXYCoord(:,:,2) = norma.*sin(alpha-theta_ref);

x_coord_circ = [];
y_coord_circ = [];

% We loop over all the buildings of the database
for n_building = 1:NumberBuildings
    % We obtain the minimum and maximum x and y coordinates of each
    % building to circumscribe each inside a rectangle
    x_max(n_building) = max(BuildingXYCoord(n_building,BuildingXYCoord(n_building,:,1)~=0,1));
    x_min(n_building) = min(BuildingXYCoord(n_building,BuildingXYCoord(n_building,:,1)~=0,1));
    y_max(n_building) = max(BuildingXYCoord(n_building,BuildingXYCoord(n_building,:,2)~=0,2));
    y_min(n_building) = min(BuildingXYCoord(n_building,BuildingXYCoord(n_building,:,2)~=0,2));
    x_coord_circ = [x_coord_circ;
                    x_min(n_building) x_min(n_building);
                    x_max(n_building) x_max(n_building);
                    x_min(n_building) x_max(n_building);
                    x_min(n_building) x_max(n_building)];
    y_coord_circ = [y_coord_circ;
                    y_min(n_building) y_max(n_building);
                    y_min(n_building) y_max(n_building);
                    y_min(n_building) y_min(n_building);
                    y_max(n_building) y_max(n_building)];
end

% for n_segment = 1:NumberBuildings*4
%     plot(x_coord_circ(n_segment,:),y_coord_circ(n_segment,:),'r','LineWidth',2)
% end

% Once all the segments have been rotated accordingly, the mean average and
% widths of the blockages are obtained
% We obtain the angles of each segment with the current coordinates system
theta = angle(exp(1i*(angle(x_coord_circ(:,2)+1i*y_coord_circ(:,2)-(x_coord_circ(:,1)+1i*y_coord_circ(:,1))))));

% Finally, we obtain the averages of the lengths and widths of the lines
% and their densities
L = sqrt((x_coord_circ(:,1)-x_coord_circ(:,2)).^2+(y_coord_circ(:,1)-y_coord_circ(:,2)).^2);
E_L = mean(L);

% Obtain the CDF
[F_L_emp,l_emp_cdf] = ecdf(L);

% We rotate the resulting vectors.
F_L_emp = F_L_emp.';
l_emp_cdf = l_emp_cdf.';

% We finally obtain the densities of segments in the scenario
A_T = (HighX-LowX)*(HighY-LowY);

N_segments = length(theta);
lambda_b = N_segments/A_T;

end

