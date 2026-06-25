fn factorial(n:usize)->usize{
(1..=n).product()
}

fn n_p_r(n:usize,r:usize)->usize{
factorial(n)/factorial(n-r)
}

fn n_c_r(n:usize,r:usize)->usize{
n_p_r(n,r)/factorial(r)
}

