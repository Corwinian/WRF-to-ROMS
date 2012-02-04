/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package gribtonetcdfconvertor;

import LRFD.Common.GeoRectangle;
import LRFD.Db.NetCDFFile.NetCDFOperator;
import java.io.File;
import java.io.IOException;
import java.util.*;
import ucar.ma2.DataType;
import ucar.ma2.Index;
import ucar.nc2.*;

/**
 *
 * @author corwin
 */
public class RomsTopLavel 
{
	
	
	public class RomsGrid
	{
		NetcdfFile existGrid;
		
		private String lon_rho_name =  "lon_rho",
				lat_rho_name = "lat_rho",
				lon_psi_name = "lon_psi",
				lat_psi_name = "lat_psi",
				lon_u_name = "lon_u",
				lat_u_name = "lat_u",
				lon_v_name = "lon_v",
				lat_v_name = "lat_v";
		
		public class grid
		{
			private double[][] lat;
			private double[][] lon;
			
			public String lon_name, lat_name;
			
			private boolean mCurvilinear = false;
			
			private double[][] loadCoord(String coordinate, boolean isLon) throws IOException
			{
				//NetCDFOperator.Read2DFieldFromFile(inData[t], "Temperature_surface")
				
				Variable var=existGrid.findVariable(coordinate);
				List<Dimension> dims=var.getDimensions();
				
				if (dims.size() != 2)				
					throw new IOException("unexpected lat dimension!");
				
				ucar.ma2.Array varar=var.read();
				
				Index ind=varar.getIndex();
				ind.set(0, 0);
				
				int jn = dims.get(isLon ? 0 : 1 ).getLength();
				int kn = dims.get(isLon ? 1 : 0 ).getLength();
				
				double res[][] = new double[isLon ? jn : kn][isLon ? kn : jn];
				
				for (int j=0; j < jn; j++)
				{
					for (int k=0; k < kn; k++)
					{
						res[isLon ? j : k][isLon ? k : j] =varar.getDouble(ind);
						ind.incr();//set(i,j,k);
						//System.out.print(String.format("%.2f ", SSTar.getDouble(ind)));
					}
				}
			
				return res;
			}
			
			public double []getLat(){ return lat[0];}
			public double []getLon(){ return lon[0];}
			
			public grid (String lonName, String latName) throws IOException
			{
				lon = loadCoord(lon_name = lonName, true);
				lat = loadCoord(lat_name = latName, false);
			}
			
			public  Map<String, Dimension> createDimension(NetcdfFileWriteable cdf) throws IOException
			{
				HashMap<String, Dimension> res = new HashMap();
		
				Dimension lonDim=cdf.addDimension(lon_name, lon[0].length),
						latDim=cdf.addDimension(lat_name, lat[0].length);
		
				Dimension[] lodim=new Dimension[1];
				lodim[0]=lonDim;
				cdf.addVariable(lon_name, DataType.DOUBLE, lodim);
				
				Dimension[] latim=new Dimension[1];
				latim[0] = latDim;
				cdf.addVariable(lat_name, DataType.DOUBLE, latim);
				
				res.put(lat_name, latDim);
				res.put(lon_name, lonDim);
				return res;
			}
		};
		
		public  RomsGrid.grid rho, psi, u, v;
		
		RomsGrid(String road) throws IOException
		{
			existGrid = NetcdfFile.open(road);
			
			rho = new RomsGrid.grid( "lon_rho", "lat_rho" );
			psi = new RomsGrid.grid("lon_psi", "lat_psi" );
			u = new RomsGrid.grid("lon_u", "lat_u" );
			v = new RomsGrid.grid("lon_v", "lat_v" );
		}
		
		public GeoRectangle getRectangle() throws Exception
		{
			return new GeoRectangle( (float)u.getLat()[0] - 4, 
					(float)u.getLat()[u.getLat().length-1] + 4,
					(float)u.getLon()[0] - 4,
					(float)u.getLon()[u.getLon().length-1] + 4);
		}
		
		public Map<String, Dimension> createDimension(NetcdfFileWriteable cdf) throws IOException
		{
			HashMap<String, Dimension> res = new HashMap();
			
			res.putAll(rho.createDimension(cdf));
			res.putAll(psi.createDimension(cdf));
			res.putAll(u.createDimension(cdf));
			res.putAll(v.createDimension(cdf));
			
			return res;
		}
		
		public double[] getParam(String name)
		{
			switch(name)
			{
				case "lon_rho":
					return rho.getLon();
				case "lat_rho":
					return rho.getLat();
				case "lon_psi":
					return psi.getLon();
				case "lat_psi":
					return psi.getLat();
				case "lon_u":
					return u.getLon();
				case "lat_u":
					return u.getLat();
				case "lon_v":
					return v.getLon();
				case "lat_v":
					return v.getLat();
				default:
					return null;
			}
		}
	};
	
	class RomsVariable
	{
		public String name;
		public String timeName, lonName, latName;
		
		public RomsVariable(String Name, String time, String type)
		{
			name = Name;
			timeName = time;
			lonName = String.format("lon_%s",type);
			latName = String.format("lat_%s",type);
		}
	};
	
	static NetcdfFileWriteable cdf;
	static String dstFile;
	
	RomsGrid grid;
	
	Map<Integer, RomsVariable> resVals;
	
	public RomsTopLavel(String dst_file, String gridFile) throws IOException
	{
		grid = new RomsGrid(gridFile);
		resVals = new HashMap<Integer, RomsVariable>(24);
		dstFile = dst_file;
		
		resVals.put(1, new RomsVariable("Temperature", "time", "u"));
		resVals.put(7, new RomsVariable("Geopotential_height", "time", "u"));
		
		resVals.put(13, new RomsVariable("Potential_temperature", "time", "u"));
		resVals.put(51, new RomsVariable("Specific_humidity", "time", "u"));
		resVals.put(52, new RomsVariable("Relative_humidity", "time", "u"));
		
		resVals.put(11, new RomsVariable("SST", "time", "u"));//Temperature surfase
		
		resVals.put(33, new RomsVariable("u_wind", "time", "u"));
		resVals.put(34, new RomsVariable("v_wind", "time", "u"));
		
		resVals.put(204, new RomsVariable("Downward_short_wave_flux", "time", "u"));
//		resVals.put(205, new RomsVariable("Downward_long_wave_flux", "time", "u"));
//		
//		resVals.put(211, new RomsVariable("Upward_short_wave_flux", "time", "u"));
//		resVals.put(212, new RomsVariable("Upward_long_wave_flux", "time", "u"));
		
		resVals.put(124, new RomsVariable("Zonal_momentum_flux", "time", "u"));
		resVals.put(125, new RomsVariable("Meridional_momentum_flux", "time", "u"));
		resVals.put(122, new RomsVariable("Sensible_heat_flux", "time", "u"));
		
		resVals.put(155, new RomsVariable("Ground_heat_flux", "time", "u"));
		resVals.put(121, new RomsVariable("Latent_heat_flux", "time", "u"));
		//resVals.put(172, new RomsVariable("Downward_long_wave_flux", "time", "u"));
		
		resVals.put(92, new RomsVariable("Ice_thickness", "time", "u"));
		resVals.put(91, new RomsVariable("Ice_concentration_ice1no_ice0", "time", "u"));
		
		resVals.put(154, new RomsVariable("Land_Surface_Precipitation_Accumulation_LSPA", "time", "u"));
		resVals.put(81, new RomsVariable("Land_cover_land1sea0", "time", "u"));
		resVals.put(57, new RomsVariable("Evaporation", "time", "u"));
		
		
		createFile();
}
	
	public RomsGrid.grid  getGridForVariable(int varNum)
	{
		return grid.u;
	}
	
	public void writeField(int fieldNum, double[][][] data)
	{
		NetCDFOperator.writeFieldToNetCDF(dstFile, resVals.get(fieldNum).name, data);
	}
	
	public void createFile() throws IOException
	{
		try
		{
			double [] time = {3.0};
			cdf = NetcdfFileWriteable.createNew(dstFile);
			
			Dimension timeDim = cdf.addDimension("time", time.length);
			Map<String,Dimension> Dimensions = grid.createDimension(cdf);
			
			Attribute att = new Attribute("missing_value", -9999);
			
			for(Iterator<Integer> i = resVals.keySet().iterator(); i.hasNext();)
			{
				int num = i.next();
				RomsVariable var = resVals.get(num);
				Dimension[] todim = new Dimension[3];
				
				todim[0] = timeDim;				
				todim[1] = Dimensions.get(var.latName); 
				todim[2] = Dimensions.get(var.lonName); 
				
				Variable varr = cdf.addVariable(var.name, DataType.DOUBLE, todim);
				
				varr.addAttribute(new Attribute("GRIB_param_number", num));
				varr.addAttribute(new Attribute("GRIB_param_name", var.name));
				
				varr.addAttribute(att);
			}
			
			cdf.create();
			cdf.close();
			
			for(Iterator<String> i = Dimensions.keySet().iterator(); i.hasNext(); i.next())
			{
				NetCDFOperator.writeTimeToNetCDF(dstFile, i.toString(), grid.getParam(i.toString()));
			}
			
			NetCDFOperator.writeTimeToNetCDF(dstFile, "time", time);
		}
		catch (Exception er)
		{
			er.printStackTrace();
		}
	}
}
