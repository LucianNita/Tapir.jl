function get_T_guess(old_mesh,new_N)
    T_new=zeros(new_N+1);
    old_N=length(old_mesh)-1;
    interval_lengths=diff(old_mesh);
    add_for_each=floor(Int,(new_N-old_N)/old_N);
    left_to_add=new_N-add_for_each*old_N;
    if left_to_add>=old_N
        error("Division and truncation went wrong")
    end
    return T_new;
end

function incremented_T_guess(old_mesh,new_N)
    T_new=zeros(new_N+1);
    old_N=length(old_mesh)-1;
    if old_N+1!=new_N
        warning("step larger than 1, redirecting to the main function")
        return get_T_guess(old_mesh,new_N)
    end
    interval_lengths=diff(old_mesh);
    ~, index = findmax(interval_lengths);
    T_new[1:index]=copy(old_mesh[1:index]); #deep?
    T_new[index+1]=0.5*(old_mesh[index]+old_mesh[index+1]);
    T_new[index+2:end]=copy(old_mesh[index+1:end]);
    return T_new;
end