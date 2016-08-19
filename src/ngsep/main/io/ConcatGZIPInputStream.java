/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.main.io;

import java.io.IOException;
import java.io.InputStream;
import java.io.PushbackInputStream;
import java.util.zip.GZIPInputStream;

/**
 * Borrowed from Trimmomatic. See http://www.usadellab.org/cms/?page=trimmomatic
 */
public class ConcatGZIPInputStream extends InputStream
{
	private PushbackInputStream source;
	private GZIPHelperInputStream gzIn;

	public ConcatGZIPInputStream(InputStream in) throws IOException
	{
		source = new PushbackInputStream(in, 1024);
		nextGzipInputStream();
	}

	private void nextGzipInputStream() throws IOException
	{
		boolean more=false;
		
		if((gzIn!=null)&&(gzIn.pushbackUnused()>0))
				more=true;
	
		if(!more)
			{
			int r=source.read();
			if(r!=-1)
				{
				source.unread(r);
				more=true;
				}
			}
		
		if(more)
			gzIn=new GZIPHelperInputStream(source);
		else
			gzIn=null;
	}
		
	@Override
	public void close() throws IOException
	{
		gzIn=null;
		source.close();
	}

	@Override
	public int read() throws IOException
	{
		int res=-1;
		
		while(res==-1 && gzIn!=null)
			{
			res=gzIn.read();
			if(res==-1)
				nextGzipInputStream();
			}
	
		/*
		if(gzIn==null)
			return -1;
	
		int res=gzIn.read();
		if(res==-1)
			{
			nextGzipInputStream();
			if(gzIn==null)
				return -1;
			else
				res=gzIn.read();
			}
	*/
		
		return res;
	}

	@Override
	public int read(byte[] b, int off, int len) throws IOException
	{
		int res=-1;
	
		while(res==-1 && gzIn!=null)
			{
			res=gzIn.read(b,off,len);
			if(res==-1)
				nextGzipInputStream();
			}
	
	/*
		if(gzIn==null)
			return -1;
	
		int res=gzIn.read(b, off, len);
		if(res==-1)
			{
			nextGzipInputStream();
			if(gzIn==null)
				return -1;
			else
				res=gzIn.read(b, off, len);
			}
		*/
	
		return res;
	}

	@Override
	public int read(byte[] b) throws IOException
	{
		int res=-1;
	
		while(res==-1 && gzIn!=null)
			{
			res=gzIn.read(b);
			if(res==-1)
				nextGzipInputStream();
			}
	
	/*
		if(gzIn==null)
			return -1;
	
		int res=gzIn.read(b);
		if(res==-1)
			{
			nextGzipInputStream();
			if(gzIn==null)
				return -1;
			else
				res=gzIn.read(b);
			}	
		*/

		return res;	
	}

	private class GZIPHelperInputStream extends GZIPInputStream
	{
		private GZIPHelperInputStream(InputStream in) throws IOException
		{
			super(in);
		}
		
		private int pushbackUnused() throws IOException
		{
			int amount=inf.getRemaining()-8;
			if(amount>0)
				source.unread(buf, len-amount, amount);
			
			return amount;
		}
	}

}
