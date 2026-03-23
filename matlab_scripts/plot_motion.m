par = load('sub-20131\ses-1\func\rp_sub-20131_ses-1_task-TAU1_run-2_bold.txt'); % Load your specific file name here
figure; % Open a new figure window
subplot(2,1,1); % Top plot for translations
plot(par(:, 1:3));
ylabel('translation (mm)');
legend({'x translation', 'y translation', 'z translation'});
subplot(2,1,2); % Bottom plot for rotations
plot(rad2deg(par(:, 4:6))); % Convert radians to degrees for better display
ylabel('rotation (degrees)');
xlabel('image (scan number)');
legend({'pitch', 'roll', 'yaw'});