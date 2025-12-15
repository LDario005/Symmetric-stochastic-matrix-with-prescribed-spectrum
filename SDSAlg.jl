module SDSAlg

    """
        return x if x is positive, 0 otherwise
    """
    MakeOh(x) = max(0.0, x)

    """
            initialize a symmetric matrix with prescribed eigenvalues
    """
    @inline function initSM(λ::Vector{Float64})
        
        # Initialize a random symmetric matrix with the given eigenvalues 
        λ=sort(λ, rev=true)
        N=length(λ)
        A = rand(0:99, N, N)
        Y = 0.5 * (A + A')
        eigvals_Y, eigvecs_Y = eigen(Y, sortby=-)
        V = eigvecs_Y
        Λ = Diagonal(λ)
        X = V * Λ * V'
        return 0.5 * (X + X'), MakeOh.(real.(X))  # Symmetrize
    end

    """
        iteration to obtain symmetric matrix with prescribed eigenvalues
    """

    @inline function projection_D( λ::Vector{Float64}, Y_new)
        eigvals_Y, eigvecs_Y = eigen(Y_new, sortby=-)
        V=eigvecs_Y
        X = V * Diagonal(λ) * inv(V)
    X = 0.5 * (X + X')  # Symmetrize again
        Y_new = MakeOh.(real.(X))
        return X, Y_new
    end

    """
        complete algorithm to obtain symmetric matrix with prescribed eigenvalues
    """
    @inline function return_SM(λ::Vector{Float64}; eps::Float64=1e-15)
        # return a random symmetric matrix with the given eigenvalues
        X, Y_new = initSM(λ)
        n=0
        while norm(X - Y_new, 2) > eps
        X, Y_new=projection_D(λ,  Y_new)
        # Y_new=Dykstra_algorithm(X)
        # Apply projection
        n+=1
        n>100_000 && break
        end
        return Y_new
    end
        """
        complete algorithm to project to find the neares stochastic symmetric matrix
    """
    @inline function Dykstra_algorithm(Y::Matrix{Float64}; eps::Float64=1e-15)
        # Initialize the symmetric matrix Y
        N= size(Y, 1)
        X=zeros(N, N)  # Initialize X as a zero matrix
        I2=zeros(N, N)
        J=ones(N, N)/N  # Initialize J as a matrix of ones
        W=I(N)-J
        n=1
        # Iterative loop until convergence
        while (norm(Y - X, 2) > eps)
            X = W * Y * W + J
            Y=MakeOh.(X-I2)
            I2=Y-X+I2
            
        end
        return MakeOh.(Y)  # Return the final matrix Y
    end

    """
        Complete procedure to find a symmetric stochastic matrix with spectrum γ;
        γ must contain the value 1 and all its entries must be real.
    """
    function Rammal_procedure(λ=Vector{Float64}; eps::Float64=1e-6)
        @assert 1.0 in λ "λ must contain the value 1"
        @assert  sum(isa.(λ, Real))==length(λ) "the spectrum must be real"
        A=@time return_SM(λ; eps=eps)
        X=zeros(size(A))
        Y=copy(A)
        n=1
        prog = ProgressUnknown(desc="Matrix construction:", spinner=true)
        while norm(X - Y, 2)≥eps
            next!(prog, spinner="🌑🌒🌓🌔🌕🌖🌗🌘")
            eigenval, V= eigen(Y, sortby=-)
            X = V * Diagonal(λ) * inv(V)
            Y=Dykstra_algorithm(X; eps=eps)
            n=n+1
            n>1e5 && break
        end
        return Y
        finish!(prog)
    end
        
end 
