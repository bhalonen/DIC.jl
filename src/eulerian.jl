function find_eulerian(result::DIC_Output{Lagrangian}, time_table::Vector{<:Real}; sample_count::Integer = 500)
    samples = map(idx->(time = time_table[idx%length(time_table)+1], point = get_roi_point(result.roi)), 1:sample_count)
    cost_function(polynomial_coeff) = cost_of_inverse_function(polynomial_coeff, result, samples)
    @time optim_result = optimize(cost_function, zeros(length(result.u_transform.a)+length(result.v_transform.a)), NelderMead())
    u_params=Optim.minimizer(optim_result)[1:length(result.u_transform.a)]
    v_params=Optim.minimizer(optim_result)[length(result.u_transform.a)+1:end]
    DIC_Output{Eulerian}(Polynomial(u_params,result.u_transform.x),Polynomial(v_params,result.v_transform.x),result.roi)
end

function back_cast(point::Point, time::Real, u_inv::Polynomial, v_inv::Polynomial)
    return (x = point.x+u_inv(point.x,point.y, time), y = point.y+v_inv(point.x,point.y, time))
end

function cost_of_inverse_function(polynomial_coeff::AbstractVector,result::DIC_Output{Lagrangian}, samples::AbstractVector{<:NamedTuple})
    u_inv = Polynomial(polynomial_coeff[1:length(result.u_transform.a)],result.u_transform.x)
    v_inv = Polynomial(polynomial_coeff[length(result.u_transform.a)+1:end],result.v_transform.x)
    sum(samples) do sample
        resampled_point = back_cast(position(sample.point.y, sample.point.x, sample.time, result), sample.time, u_inv, v_inv)
        sqrt((sample.point.x - resampled_point.x)^2 + (sample.point.y - resampled_point.y)^2)
    end
end