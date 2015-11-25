all : modshift jktebop jktld

modshift :
	g++ -Wno-unused-result -O3 -o modshift modshift-v4.cpp

jktebop :
	gfortran -O3 -o jktebop jktebop-dave.f

jktld :
	gfortran -O3 -o jktld jktld-v3.f

clean:
	rm bin-noerr modshift jktebop jktld
