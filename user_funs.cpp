#include"user_funs.h"
#include<iomanip>

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera wartoœæ funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wspó³rzêdne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera wartoœæ funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz¹tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si³y dzia³aj¹cy na wahad³o oraz czas dzia³ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi¹zujemy równanie ró¿niczkowe
	int n = get_len(Y[0]);									// d³ugoœæ rozwi¹zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad³a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// wartoœæ funkcji celu (ud1 to za³o¿one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po³o¿enia to prêdkoœæ
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z prêdkoœci to przyspieszenie
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2) {
	return -cos(0.1 * x(0)) * exp(-pow((0.1 * x(0) - 2 * 3.14), 2)) + 0.002 * pow((0.1 * x(0)), 2);
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
	matrix dY(Y);

	double a = 0.98;
	double b = 0.63;
	double g = 9.81;
	double PA = 0.5;
	double PB = 1;
	double DB = 0.00365665;
	double Fin = 0.01;
	double Tin = 20;
	double TA = 90.0;
	double DA = ud1(0);

	double VA = Y(0);
	double VB = Y(1);
	double TB = Y(2);

	double FAout = VA > 0 ? a * b * DA * sqrt(2 * g * VA / PA) : 0;
	double FBout = VB > 0 ? a * b * DB * sqrt(2 * g * VB / PB) : 0;

	dY(0) = -FAout;
	dY(1) = FAout + Fin - FBout;
	dY(2) = FAout / VB * (TA - TB) + Fin / VB * (Tin - TB);

	return dY;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2)
{
	matrix y(1, 1);
	matrix y0 = matrix(3, new double[3] { 5, 1, 20 });

	matrix* y_ptr = solve_ode(df1, 0, 1, 2000, y0, ud1, x);

	int n = get_len(y_ptr[0]);

	double max = y_ptr[1](0, 2);

	for (int i = 0; i < n; i++)
		if (y_ptr[1](i, 2) > max)
			max = y_ptr[1](i, 2);

	cout << std::setprecision(20) << "Max = " << max << endl;

	y(0) = abs(max - 50.0);

	y_ptr[0].~matrix();
	y_ptr[1].~matrix();

	return y;
}
