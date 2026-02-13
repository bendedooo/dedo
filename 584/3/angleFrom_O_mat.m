function ang_out = angleFrom_O_mat(O_mat)

    ang_out(1,:) = atan2(O_mat(2,3:3:end),O_mat(3,3:3:end));
    ang_out(2,:) = -asin(O_mat(1,3:3:end));
    ang_out(3,:) = atan2(O_mat(1,2:3:end),O_mat(1,1:3:end));
end