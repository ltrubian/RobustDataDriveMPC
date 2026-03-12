function g = calculate_gamma(P, theta)
    % gamma = 0.5 * [tr((I - theta*P)^-1 - I) + ln(det(I - theta*P))]
    n = size(P, 1);
    In = eye(n);
    M = In - theta * P;
    g = 0.5 * (trace(inv(M) - In) + log(det(M)));
end