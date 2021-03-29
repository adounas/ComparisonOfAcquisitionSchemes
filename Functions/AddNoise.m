function ksp = AddNoise(ksp,std)

ksp = ksp + std*randn(size(ksp)) + 1i*std*randn(size(ksp));

end