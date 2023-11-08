using AQFED, AQFED.Math, AQFED.Basket, Printf, DataFrames

data = DataFrame(n=Int[], strike=Float64[], sigma=Float64[], rho=Float64[], price=Float64[], Legendre=Float64[], Chebyshev1=Float64[], Chebyshev2=Float64[], TanhSinh=Float64[], Simpson=Float64[], Deelstra=Float64[])

nAsset = 2
	rhoMin = -1.0
	rhoMax = 1.0
	rhoSteps = 10
    pSimp = AQFED.Basket.QuadBasketPricer(AQFED.Math.Simpson(1024*64))
    rhos = [-0.9,-0.5,0.5,0.9]
 	for rho=rhos
		r = 0.05
        q = 0.0
        		σ = 0.2
		tte = 1.0
		strike = 100.0
        nWeights = 2
		discountFactor = exp(-r * tte)
        weights = zeros(Float64, nWeights)
        tvar = zeros(Float64, nWeights)
        forward = zeros(Float64, nWeights)
        spot=zeros(nWeights)
        for i = eachindex(weights)
            weights[i] = 1.0 / (nWeights)
            tvar[i] = σ^2 * tte
            spot[i] = 100.0
            forward[i] = spot[i] * exp((r - q) * tte)
        end
    
        correlation = [
            1.0 rho 
            rho 1.0 
        ]
        price = priceEuropean(pSimp, true, strike, discountFactor, spot, forward, tvar, weights, correlation,isSplit=true)
        for n=4:128
        pLeg = AQFED.Basket.QuadBasketPricer(AQFED.Math.GaussLegendre(n))
        pSin = AQFED.Basket.QuadBasketPricer(AQFED.Math.TanhSinh(n,1e-16))#,2*4.026/(n-1))) # 2tmax/(n-1)
        pCheb1 = AQFED.Basket.QuadBasketPricer(AQFED.Math.Chebyshev{Float64,1}(n))
        pCheb2 = AQFED.Basket.QuadBasketPricer(AQFED.Math.Chebyshev{Float64,2}(n))
        pd = AQFED.Basket.DeelstraBasketPricer(3,3,AQFED.Math.GaussLegendre(n))    
        priceLeg = priceEuropean(pLeg, true, strike, discountFactor, spot, forward, tvar, weights, correlation,isSplit=true,isLognormal=true)
        priceCheb1 = priceEuropean(pCheb1, true, strike, discountFactor, spot, forward, tvar, weights, correlation,isSplit=true,isLognormal=true)
        priceCheb2 = priceEuropean(pCheb2, true, strike, discountFactor, spot, forward, tvar, weights, correlation,isSplit=true,isLognormal=true)
        priceSin = priceEuropean(pSin, true, strike, discountFactor, spot, forward, tvar, weights, correlation,isSplit=true,isLognormal=true)
        priceSim = priceEuropean(AQFED.Basket.QuadBasketPricer(AQFED.Math.Simpson(n)), true, strike, discountFactor, spot, forward, tvar, weights, correlation,isSplit=true,isLognormal=true)
        priceD = priceEuropean(pd, true, strike, discountFactor, spot, forward, tvar, weights, correlation)
        @printf("%d %.0f %.2f %.1f %.2f %f %.2e %.2e %.2e %.2e %.2e\n", n, strike, r, σ, rho, price, priceLeg-price,priceCheb1-price, priceCheb2-price,priceSin-price,  priceD-price)
            push!(data, [n, strike,  σ, rho, price, priceLeg-price,priceCheb1-price, priceCheb2-price,priceSin-price,  priceSim-price,priceD-price])
    end
		# if math.Abs(priceRef-refValues[ir]) > 1e-8 {
		# 	t.Errorf("error too large at %d, expected %.8f was %.8f", ir, refValues[ir], priceRef)
		# }
    end
    theme(:juno)
    allPlots = []
    for rho=rhos
    df50 = data[data.rho .== rho,:]
    p1 = plot(df50.n,abs.(df50.Legendre).+1e-20,label="Legendre")
    plot!(p1,df50.n,abs.(df50.Chebyshev1).+1e-20,label="Chebyshev1")
    plot!(p1,df50.n,abs.(df50.Chebyshev2).+1e-20,label="Chebyshev2")
    plot!(p1,df50.n,abs.(df50.TanhSinh).+1e-20,label="TanhSinh")
    plot!(p1,df50.n,abs.(df50.Simpson).+1e-20,label="Simpson")
    plot!(p1,xlab="N",ylab="Absolute error",ylim=[1e-16,1.0],yscale=:log10,size=(640,480))
    push!(allPlots,p1)
    end
    plot(allPlots[1],allPlots[2],allPlots[3],allPlots[4],layout=(4,1),size=(600,1600),left_margin=12Plots.mm)
    plot!(yticks=[1,1e-2,1e-4,1e-6,1e-8,1e-10,1e-12,1e-14])
    savefig(string("~/Dev/quad_rho_all.png"))
    