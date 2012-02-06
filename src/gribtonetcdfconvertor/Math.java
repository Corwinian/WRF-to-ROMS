/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package gribtonetcdfconvertor;

import java.lang.Math;
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

	public Complex sqrt() 
	{
		return new	Complex(Math.sqrt(mod()) * Math.cos(arg()), Math.sqrt(mod()) * Math.sin(arg()));
	}
		
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

	public static Complex abs(Complex z){ return new Complex(Math.abs(z.re), Math.abs(z.im));}
	public Complex abs(){ return Complex.abs(this);}
	
	public Complex plus(double z){ return plus(this, new Complex(z));}
	public Complex plus(Complex z){ return Complex.plus(this, z);}
	public static Complex plus(double a, Complex b){ return plus (new Complex(a), b);}
	public static Complex plus(Complex a, double b){ return plus (a, new Complex(b));}	
	public static Complex plus(Complex a, Complex b)
	{ return new Complex(a.re + b.re, a.im + b.im);}
	
	public Complex minus(double z){ return minus(this, new Complex(z));}
	public Complex minus(Complex z){ return Complex.minus(this, z);}
	public static Complex minus(double a, Complex b){ return minus (new Complex(a), b);}
	public static Complex minus(Complex a, double b){ return minus (a, new Complex(b));}	
	public static Complex minus(Complex a, Complex b)
	{return  new Complex(a.re - a.re, a.im - b.im);}
	
	public Complex asterisk(double z) {return asterisk(new Complex(z));}
	public Complex asterisk(Complex z){return asterisk(this, new Complex(z));} 
	public static Complex asterisk(double a, Complex b){ return asterisk (new Complex(a), b);}
	public static Complex asterisk(Complex a, double b){ return asterisk (a, new Complex(b));}	
	public static Complex asterisk(Complex a, Complex b)
	{
		return new Complex( a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
	}
	
	public Complex slash(double z){ return slash(this, new Complex(z));}
	public Complex slash(Complex z){ return Complex.slash(this, z);}
	public static Complex slash(double a, Complex b){ return slash (new Complex(a), b);}
	public static Complex slash(Complex a, double b){ return slash (a, new Complex(b));}	
	public static Complex slash(Complex a, Complex b)
	{ 
		double m = b.mod();  
		return new Complex((a.re * b.re - a.im * b.im) / m, (a.im * b.re - a.re * b.im) / m); 
	}
};

class Vectorr
{
	
};