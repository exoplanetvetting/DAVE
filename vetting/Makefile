all : modshift jktebop jktld

modshift : modshift.cpp
	g++ -Wno-unused-result -O3 -o modshift modshift.cpp

jktebop :
	gfortran -O3 -o jktebop jktebop-dave.f

jktld :
	gfortran -O3 -o jktld jktld-v3.f

clean:
	rm bin-noerr modshift jktebop jktld
