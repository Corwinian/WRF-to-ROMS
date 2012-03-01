/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package gribtonetcdfconvertor;

import java.io.*;
import LRFD.Common.*;

import LRFD.DataStructures.*;
import LRFD.Methods.Interpolation;
import LRFD.Db.NetCDFFile.*;

import ucar.nc2.*;
import ucar.ma2.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Iterator;
import java.util.HashMap;
import javax.swing.plaf.synth.Region;
import sun.misc.Sort;
import java.lang.Math;

import java.math.MathContext;


 /**
 *
 * @author corwin
 */
public class GribToNetCDFConvertor 
{

	private static String latName = "lat";
	private static String lonName = "lon";
	private static String timeName = "time";
	private static String time1Name = "time1";
	private static String levelName = "height_above_ground";
	private static String level1Name = "height_above_ground1";
	
	private static VariablesNums[] neededValues = {
		VariablesNums.Pressure,//Surface pressure
		VariablesNums.Geopotential_height,
		//Terrain height
		VariablesNums.Potential_temperature,//Skin potential temperature
		VariablesNums.Specific_humidity,//Skin specific humidity
		VariablesNums.Relative_humidity,//Skin Relative humidity
		VariablesNums.SST,//Skin temperature
		VariablesNums.u_wind,//10 M u component wind
		VariablesNums.v_wind,//10 M v component wind
	//	VariablesNums.dsrad ,//Incoming surface shortwave radiation — time-averaged
		VariablesNums.svstr,//Surface u wind stress
		VariablesNums.sustr,//Surface v wind stress
		VariablesNums.Sensible_heat_flux,//Surface sensible heat flux — time-averaged
		VariablesNums.Ground_heat_flux,//Ground heat flux — time-averaged
		VariablesNums.Latent_heat_flux,//Surface latent heat flux — time-averaged
	//	172,//
		VariablesNums.Ice_concentration_ice1no_ice0,//Sea ice mask
		VariablesNums.Ice_thickness,//Ice thickness
		VariablesNums.Land_cover_land1sea0,//Land sea mask (land=1, sea=0)
		VariablesNums.Land_Surface_Precipitation_Accumulation_LSPA,//Accumulated land surface model precipitation
		VariablesNums.Evaporation,//Accumulated surface evaporation
	};
	
	private static Map<VariablesNums, String> getVariablesByNums(NetcdfFile cdf) throws IOException
	{
		List<Variable> variables = cdf.getVariables();
		
		Map<VariablesNums, String> map = new HashMap<>();
		/*я не нашел как сделать через итераторы, поэтому сделал так*/
		for(int i = 0; i < variables.size(); ++i)
		{ 
			Attribute atr = variables.get(i).findAttribute("GRIB_param_number");
			
			try
			{
				//FIXME:потом исправить
				if (atr != null)
				{
					if ((Integer)atr.getValue(0) == 11 &&
						variables.get(i).getName().equals("Temperature"))
						map.put(VariablesNums.Temperature, variables.get(i).getName()); //сделал из-за то что у температуры и температуры поверхности одинаковые номера
					else			
						map.put(VariablesNums.getType((Integer)atr.getValue(0)), variables.get(i).getName());
				}
			}
			catch(Exception ex)
				{}
		}
		
		return map;
	}
	
	private static double[] loadCoords(NetcdfFile cdf, String coordinate) throws IOException
	{
		Variable latitude = cdf.findVariable(coordinate);
		ucar.ma2.Array lats = latitude.read();
		int n = (int)lats.getSize();
		double lat[] = new double[n];
		
		for (int i = 0; i < n; i++)
		{
			lat[i]=lats.getDouble(i);
		}
		return lat;
	}
	
	private static double[] loadCoords(NetcdfFile cdf, String coordinate, double min, double max) throws IOException, Exception
	{
		return loadCoords(cdf, coordinate, min, max, false);
	}
	
	private static double[] loadCoords(NetcdfFile cdf, String coordinate, double min, double max, boolean inv) throws IOException, Exception
	{
		Variable var = cdf.findVariable(coordinate);
		
		Range range=getRange(var, min, max, false);
		return selectRange(var, range, inv);
	}
	
	public static Range getRange(Variable Var, double minValue, double maxValue) throws Exception
	{
		return getRange(Var, minValue, maxValue, false);
	}
	
	public static Range getRange(Variable Var, double minValue, double maxValue, boolean wide) throws Exception
	{
		ucar.ma2.Array tvar=Var.read();
		int iMin=0, iMax=0;
		if (Var.getRank()!=1) throw new Exception("Error in dimensions. Rank != 1");

		int n=Var.getDimension(0).getLength();
		double minVal=Double.MAX_VALUE;
		double maxVal=-Double.MAX_VALUE;

		for (int i=0; i < n; i++)
		{
			double timevar = tvar.getDouble(i);
			
			if (timevar >= minValue && timevar <= maxValue)
			{
				if (timevar < minVal)
				{
					iMin=i;
					minVal=timevar;
				}
				if (timevar > maxVal)
				{
					iMax=i;
					maxVal=timevar;
				}
			}
		}
		if (iMax<iMin)
		{
			int temp=iMax;
			iMax=iMin;
			iMin=temp;
		}
		if (minVal!=minValue || maxVal!=maxValue)
		{
			iMin--;
			iMax++;
			if (iMax>=n) 
				iMax=n;
			if (iMin<0) 
				iMin=0;
		}
		if (wide)
		{
			iMin-=3;
			iMax+=3;
			if (iMax>=n) 
				iMax=n;
			if (iMin<0) 
				iMin=0;
		}
		return new Range(iMin, iMax);
	}
	
	public static double[] selectRange(Variable Var, Range range) throws IOException
	{ return selectRange(Var, range, false);}
	
	public static double[] selectRange(Variable Var, Range range, boolean inv) throws IOException
	{
		double res[]=new double[range.last()-range.first()+1];
		ucar.ma2.Array tvar=Var.read();
		int first=range.first();
		
		for (int i=0; i<res.length; i++)
		{
			res[i] = inv ? tvar.getDouble(res.length - 1 - i + first) : tvar.getDouble(i+first);
		}
		
		return res;
	}
	
	public static Data3DField getFieldFromSRCFile(NetcdfFile src, String FieldName, GeoRectangle boundaries, double timeMin, double timeMax) throws Exception
	{
		return getFieldFromSRCFile(src, FieldName, boundaries, timeMin, timeMax, false);
	}
	public static Data3DField getFieldFromSRCFile(NetcdfFile src, String FieldName, GeoRectangle boundaries, double timeMin, double timeMax, boolean wide) throws Exception
	{
		Variable var=src.findVariable(FieldName);
		if (!(var.getRank()==3 || (var.getRank()==4 && var.getDimension(1).getLength()==1))) 
			throw new Exception("Error in Dimensions, Rank of dimensions for variable "+var.getName()+" must be 3.");
		
		String lvName = (FieldName.equals("Specific_humidity") || FieldName.equals("Relative_humidity")	|| FieldName.equals("Temperature")) ?
				level1Name : levelName ;
		String tmName = (FieldName.equals("Evaporation") || FieldName.equals("Land_Surface_Precipitation_Accumulation_LSPA")) ?
				time1Name : timeName;
		
		int nv=var.getRank();
		Variable timeVar=src.findVariable(tmName);
		Variable lonVar=src.findVariable(lonName);
		Variable latVar=src.findVariable(latName);

		Range timeRange=getRange(timeVar, timeMin, timeMax);
		double time[]=selectRange(timeVar, timeRange);

		Range lonRange=getRange(lonVar, boundaries.north_west.lon, boundaries.south_east.lon, wide);
		double lon[]=selectRange(lonVar, lonRange);

		Range latRange=getRange(latVar, boundaries.south_east.lat, boundaries.north_west.lat, wide);
		double lat[]=selectRange(latVar, latRange);

		ArrayList<Range> lr=new ArrayList<Range>();
		lr.add(var.findDimensionIndex(tmName), timeRange);
		if (nv==4) 
		{
			Range levelRange=new Range(0, 0);
		
			lr.add(var.findDimensionIndex(lvName), levelRange);
		}
		
		lr.add(var.findDimensionIndex(latName), latRange);
		lr.add(var.findDimensionIndex(lonName), lonRange);

		ucar.ma2.Array tvar=var.read(lr);
		long nt=tvar.getSize();
		if (nt!=lon.length*lat.length*time.length)
		{
			System.out.println("Error in dimensions");
		}

		int[] dims=new int[nv];
		dims[var.findDimensionIndex(tmName)]=time.length;
		if (nv==4) 
			dims[var.findDimensionIndex(lvName)]=1;
		dims[var.findDimensionIndex(latName)]=lat.length;
		dims[var.findDimensionIndex(lonName)]=lon.length;
		
		double res[][][] = (nv != 4) ? new double[dims[0]][dims[1]][dims[2]]:new double[dims[0]][dims[2]][dims[3]];
		double scaleFactor=1;
		double addOffset=0;
		boolean kridControl=false;
		double minV=0;
		double maxV=0;
		double mis=Double.MAX_VALUE;
		double fill=Double.MAX_VALUE;
		try
		{
			Attribute SF=var.findAttribute("scale_factor");
			Attribute AO=var.findAttribute("add_offset");

			scaleFactor=SF.getNumericValue().doubleValue();
			addOffset=AO.getNumericValue().doubleValue();
			Attribute VR=var.findAttribute("valid_range");
			minV=VR.getNumericValue(0).doubleValue();
			maxV=VR.getNumericValue(1).doubleValue();
			minV=minV*scaleFactor+addOffset;
			maxV=maxV*scaleFactor+addOffset;
			kridControl=true;
			Attribute MIS=var.findAttribute("missing_value");
			Attribute FILL=var.findAttribute("_FillValue");
			mis=MIS.getNumericValue().doubleValue();
			fill=FILL.getNumericValue().doubleValue();
			//mis=mis*scaleFactor+addOffset;
			//fill=fill*scaleFactor+addOffset;
		}
		catch (Exception er){}
		
		Index ima=tvar.getIndex();
		int first=dims[0];
		int second = (nv==4) ? dims[2]:dims[1];
		int third = (nv==4) ? dims[3]:dims[2];
		
		for (int i=0; i<first; i++)
		{
			for (int j=0; j<second; j++)
			{
				for (int k=0; k<third; k++)
				{
					if (nv==3) 
						ima.set(i, j, k);
					else if (nv==4) 
						ima.set(i, 0, j, k);
					res[i][j][k]=tvar.getDouble(ima);
					
					if (res[i][j][k]==mis)
					{
						System.out.println("Missing value = "+res[i][j][k]+" in position "+i+"; "+j+"; "+k+"; replaced by "+fill);
						res[i][j][k]=fill;
					}
					
					res[i][j][k]=res[i][j][k]*scaleFactor+addOffset;
					if (kridControl)
					{					  
						if (res[i][j][k]<minV || res[i][j][k]>maxV)
						{
							System.out.println("Illegar value = "+res[i][j][k]+" in position "+i+"; "+j+"; "+k+";");
							//throw new Exception("Illegar value in position "+i+"; "+j+"; "+k+";");
						}
					}
				}
			}
		}
		return new Data3DField(res, lat, lon, time);
	}
	
	public static Data3DField InterpolateField(Data3DField field, RomsTopLavel.RomsGrid.grid grid) throws Exception
	{
		Data3DField resflx=new Data3DField();
		resflx.lat=grid.getLat();
		resflx.lon=grid.getLon();
		resflx.time=field.time;
		resflx=Interpolation.BilinearInterpolation(field, resflx);
		return resflx;
	}
	
	
	private static ArrayList<Integer> findZeros(double[] u, double[] v)
	{
		ArrayList<Integer> res = new ArrayList<Integer>();
		for (int i=0; i < u.length; ++i)
		{
			if((Math.abs(u[i]) + Math.abs(v[i])) ==0)
				res.add(i);
		}
		
		return res;
	}
	
	private static 	ArrayList<Integer> findGood(double[] u, double[] v)
	{
		ArrayList<Integer> res = new ArrayList<Integer>();
		for (int i=0; i < u.length; ++i)
		{
			if((Math.abs(u[i]) + Math.abs(v[i])) > 0)
				res.add(i);
		}
		
		return res;
	}
	
	public static double []createTemp(ArrayList<Integer> needs, double[] all)
	{
		double [] res = new double[needs.size()];
		
		for(int i =0; i < needs.size(); i++ )
		{
			res[i] = all[needs.get(i)];
		}
		return res;
	}
	
	
	private static boolean mayCalc(Complex [] c, Complex [] d1)
	{
		double max = 0;
		for (int i=0; i < d1.length; i++)
		{
			double tem = Complex.abs(c[i].minus(d1[i]));
			if (tem > max)
				max = tem;
		}
		//TODO: подумать а правильно ли проверяю, не нужноли еще проверять мнмую часть
		return max > 0.01;
	}
	
	public static ArrayList<Integer>find_c(Complex []c)
	{
		ArrayList<Integer> res = new ArrayList<>();
		for (int i=0; i < c.length; i++)
		{
			if(c[i].getRe() > 11)
				res.add(i);
		}
		//TODO: подумать а правильно ли проверяю, не нужноли еще проверять мнмую часть
		return  res;
	}
	
	private static double[][] getWstress(double []u, double[] v)
	{
		double z0 = 10, k=0.41, rho=1.25e-3, a1=1/k*Math.log(z0/10.);
		
		double []u10 = new double[u.length], v10= new double[u.length], 
				taux= new double[u.length], tauy = new double[u.length];
		
		ArrayList<Integer> izeros = findZeros(u, v);
		ArrayList<Integer> igood = findGood(u, v);
		
		u = createTemp(igood, u);
		v = createTemp(igood, v);
		
		Complex[] w = new Complex[v.length];
		double[] v0 = new double[v.length];
		
		for (int i =0; i < w.length; ++i)
		{
			w[i]= new Complex(u[i], v[i]);
			v0[i] = w[i].abs();
		}
		
		Complex[] c = new Complex[u.length], d1=new Complex[u.length], cd = new Complex[u.length];
		
		for (int i=0; i < d1.length;i++)
		{
			d1[i] = new Complex(1.0e+35);
			c[i] = new Complex(0);
		}
		
		while(mayCalc(c, d1))
		{
			c = d1.clone();
			
			for (int i=0; i < cd.length; cd[i++] = new Complex(1.205e-3));
			
			ArrayList<Integer> ind =find_c(c);
			
			if (!ind.isEmpty())
			{
				for(Iterator<Integer> i = ind.iterator(); i.hasNext();)
				{
					Integer t = i.next();
					cd[t] = Complex.asterisk(Complex.plus(0.49, Complex.asterisk(0.065, c[t])), 1.e-3);
				}
			}
			//d1=v0./(1+a1.*sqrt(cd));
			
			for (int i =0; i < d1.length; ++i)
			{
				d1[i] = Complex.slash(v0[i], Complex.plus(1, Complex.asterisk(a1, cd[i].sqrt())));
			}
		}
		
		Complex[] t = new Complex[u.length];
		
		for (int i=0; i < t.length; i++)
		{
			t[i] = Complex.asterisk(rho, cd[i].asterisk(Complex.asterisk(1.e+4, d1[i])));
		}
		
		Complex[] w10 = new Complex[u.length];
		for (int i=0; i < w10.length; i++)
		{
			w10[i] = Complex.asterisk(d1[i].slash(v0[i]), w[i]);
		}
		
		for (int i=0; i < igood.size(); i++)
		{
			u10[igood.get(i)] = w10[i].getRe();
			v10[igood.get(i)] = w10[i].getim();
		}
		
		for (int i=0; i < izeros.size(); i++)
		{
			v10[izeros.get(i)] = u10[izeros.get(i)] = 0;
		}
		
		return new double[][]{u10, v10};
	}
	
	public static Data3DField [] CalcWstress(Data3DField u, Data3DField v)
	{
		for (int i=0; i < u.data.length; ++i)
		{
			for (int c=0; c < u.data[0].length; ++c)
			{
				double [][] t = getWstress(u.data[i][c], v.data[i][c]);
				u.data[i][c] = t[0].clone();
				v.data[i][c] = t[1].clone();
			}
		}
		
		return new Data3DField[]{u, v};
	}
	
	private static Data3DField CalcShflux(Data3DField dswr, Data3DField dlwr, Data3DField uswr, Data3DField ulwr)
	{
		for (int i =0; i < dswr.data.length; ++i)
		{
			for (int c =0; c < dswr.data[i].length; ++c)
			{
				for (int j =0; j < dswr.data[i][c].length; ++j)
				{
					dswr.data[i][c][j] += dlwr.data[i][c][j] - uswr.data[i][c][j] - ulwr.data[i][c][j];
				}
			}
		}
		
		return dswr;
	}
	
	//public static void GribToNetCDFConvert(File fileIn, String gridFile, String outFile) throws IOException, Exception //, String dstFileName)
	public static void GribToNetCDFConvert(List<String>filesIn, String gridFile, String outFile) throws IOException, Exception //, String dstFileName)
	{   
		NetcdfFile cdf = null;
		try
		{ 
			double[] time = {3, 4, 15, 204, 345, 350} ;//new double[filesIn.size()];
			
//			for (int i=0; i < filesIn.size() ; ++i)
//			{
//				cdf = NetcdfFile.open(filesIn.get(i));
//				time[i] = loadCoords(cdf, timeName)[0];
//			}
			
			RomsTopLavel dest = new RomsTopLavel(outFile, gridFile, time);
			GeoRectangle gr = dest.grid.getRectangle();
			
			Map<VariablesNums, Data3DField> fields = new HashMap<>(dest.getResVals().size());
			{
				cdf = NetcdfFile.open(filesIn.get(0));
								
				Map<VariablesNums, String> variables = getVariablesByNums(cdf);

				String field = variables.get(neededValues[0]);

				Data3DField SST = getFieldFromSRCFile(cdf, field,gr, time[0], time[0]);
				SST.data = new double[time.length][SST.data[0].length][SST.data[0][0].length];
				SST.time = time.clone();
				
				for(Iterator<VariablesNums> i = variables.keySet().iterator(); i.hasNext();)
				{
					VariablesNums var = i.next();
					fields.put(var, SST.copy());
				}
			}
			
			Map<VariablesNums, String> variables = null;

			for (int c = 0; c < filesIn.size(); ++c)
			{	
				cdf = NetcdfFile.open(filesIn.get(c));
				
				variables = getVariablesByNums(cdf);

				for (Iterator<VariablesNums> i = variables.keySet().iterator(); i.hasNext();)
				{
					VariablesNums var = i.next();
					Data3DField SST = getFieldFromSRCFile(cdf, variables.get(var), gr, time[0], time[0]);
					//Data3DField SST = getFieldFromSRCFile(cdf, variables.get(var), gr, time[c], time[c]);
					fields.get(var).data[c] = SST.data[0].clone();
				}
			}

			Map<VariablesNums, Data3DField> loaded = new HashMap<>(dest.getResVals().size());
			
			{
				for(int i = 0; i < neededValues.length; ++i)
				{
					String field = variables.get(neededValues[i]);
					if (field == null)
					{
						System.out.println(String.format("Can't find variavle with number: %d", neededValues[i]));
						continue;
					}
					System.out.println(field);
					loaded.put(neededValues[i], fields.get(neededValues[i]));
				}
				
				System.out.println("shflux");
				Data3DField teemp = CalcShflux(fields.get(VariablesNums.dsrad), fields.get(VariablesNums.dlrad), 
						fields.get(VariablesNums.usrad), fields.get(VariablesNums.ulrad));
				loaded.put(VariablesNums.shflux, teemp
						);
						
				System.out.println("dQdSST");
				//WARNING: не уверен насколько ето верно что что б не дублировать код 
				//с другого модуля имеет смысл сделать вот такое вот
				CreateForcing forsing = new CreateForcing(filesIn.get(0));
				//WARNING: я очень не уверен что подставил именно те переменнные которые нужны
				loaded.put(VariablesNums.dQdSST, forsing.getdQdSST(
							fields.get(VariablesNums.Temperature),
							fields.get(VariablesNums.Specific_humidity),
							fields.get(VariablesNums.SST),
							fields.get(VariablesNums.u_wind),
							fields.get(VariablesNums.v_wind),
							forsing.getAirDensity(
								fields.get(VariablesNums.Temperature),
								fields.get(VariablesNums.Specific_humidity),
								fields.get(VariablesNums.Pressure)
							)
						));
				
				System.out.println("windStress");
				Data3DField[] windStress = CalcWstress(fields.get(VariablesNums.u_wind), fields.get(VariablesNums.v_wind));
							
				loaded.put(VariablesNums.sustr, windStress[0]);
				loaded.put(VariablesNums.svstr, windStress[1]);
			}
			
			for(Iterator<VariablesNums> i = loaded.keySet().iterator(); i.hasNext();)
			{
				VariablesNums var = i.next();
				
				Data3DField field = loaded.get(var);
				field.InverseLatIfNecessary();
				dest.writeField(var, InterpolateField(field, dest.getGridForVariable(var)).data);
			}
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
		finally
		{
			cdf.close();
		}
	}

	
   /**
	* @param args the command line arguments
	*/
	public static void main(String[] args) throws IOException, Exception
	{

//		//String fileIn = "/home/corwin/Dropbox/Учеба/Курсовик/wrfprs.000.grb";
//		String fileIn ="/home/corwin/Dropbox/Учеба/Курсовик/wrf/wrfprs_for_roms.003.grb";
//	//	String fileIn ="/home/corwin/Dropbox/Учеба/Курсовик/wrf/wrfprs_for_roms.003.nc";
//		String fileGrid = "/home/corwin/Dropbox/Учеба/Курсовик/wrf/roms_grd.nc";
//		String fileOut = "/home/corwin/Dropbox/Учеба/Курсовик/wrf/wrfprs_for_roms.003.nc";
	
// -g roms_grd.nc -i wrfprs_for_roms.003.grb  -o wrfprs_for_roms.003.nc
		List<String>filesIn =  new ArrayList<>();
		String fileIn = "";
		String fileGrid = "";
		String fileOut = "";
		
		for (int i=0; i < args.length; i+=2)
		{
			if (args[i].equals("-g") || args[i].equals("--grid"))
				fileGrid = args[i+1];
			else if (args[i].equals("-i") || args[i].equals("--input"))
				//fileIn = args[i+1];
				filesIn.add(args[i+1]);
			else if (args[i].equals("-o") || args[i].equals("--output"))
				fileOut = args[i+1];
		}
		
		File fIn = new File(fileIn);
		
		//GribToNetCDFExtractor.printMetaData(fIn);
		GribToNetCDFConvert(filesIn, fileGrid, fileOut);
	}
}
