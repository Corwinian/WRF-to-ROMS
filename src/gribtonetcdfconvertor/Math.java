/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package gribtonetcdfconvertor;

/**
 *
 * @author corwin
 */
//скопировал класс комплекс из книги, тк как я понял в java нет класса комплексных чисел О_О
//
class Complex  
{
	private static final double EPS = 1e-12; // Точность вычислений  
	private double re, im;					// Действительная и мнимая часть 
	// Четыре конструктора  

	Complex(double re, double im) { this.re = re; this.im = im; } 
	Complex(double re){ this(re, 0.0); } 
	Complex(Complex z){this(z.re, z.im) ; } 
	Complex(){ this(0.0, 0.0); } 

	// Методы доступа  
	//
	public double getRe(){return re;}  
	public double getim(){return im;}  
	public Complex getZ(){return new Complex(re, im);}  

	public void setRe(double re){this.re = re;}  
	public void setlm(double im){this.im = im;}  
	public void setZ(Complex z){re = z.re; im = z.im;} 
	
	//
	//	// Модуль и аргумент комплексного числа 
	//
	public double mod(){return Math.sqrt(re * re + im * im);}  
	public double arg(){return Math.atan2(re, im);}

	public boolean isReal(){return Math.abs(im) < EPS;} 

	// Вывод на экран 
	public void pr()
	{	System.out.println(re + (im < 0.0 ? "" : "+") + im + "i");} 

	// Переопределение методов класса Object 
	//
	public boolean equals(Complex z)
	{	return Math.abs(re - z.re) < EPS && Math.abs(im - z.im) < EPS;} 
	public String toString(){ return "Complex: " + re + " " + im;} 

	// Методы, реализующие операции +=, -=, *=, /=  
	public void add(Complex z){re += z.re; im += z.im;}  
	public void sub(Complex z){re -= z.re; im -= z.im;}  
	public void mul(Complex z)
	{
		double t = re * z.re - im * z.im;  
		im = re * z.im + im * z.re;  
		re = t; 
	} 

	public void div(Complex z)
	{ 
		double m = z.mod(); 
		double t = re * z.re - im * z.im; 
		im = (im * z.re - re * z.im) / m; 
		re = t / m;  
	} 

	// Методы, реализующие операции +, -, *, /  

	public Complex plus(Complex z){ return new Complex(re + z.re, im + z.im);}
	public Complex minus(Complex z) {return new Complex(re - z.re, im - z.im); }
	
	public Complex asterisk(Complex z)
	{	
		return new Complex( re * z.re - im * z.im, re * z.im + im * z.re);
	} 
	
	public Complex slash(Complex z)
	{  
		double m = z.mod();  
		return new Complex((re * z.re - im * z.im) / m, (im * z.re - re * z.im) / m); 
	} 
};