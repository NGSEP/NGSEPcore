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
package ngsep.main;

import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

/**
 * 
 * @author Jorge Gomez
 *
 */
public class ThreadPoolManager {
	private static final int TIMEOUT_SECONDS = 30;
	
	private int maxTaskCount;
	private int secondsPerTask=1;
	private final int numThreads;
	private ThreadPoolExecutor pool;
	
	public ThreadPoolManager(int numberOfThreads, int maxTaskCount) {
		this.pool = new ThreadPoolExecutor(numberOfThreads, numberOfThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		this.maxTaskCount = maxTaskCount;
		this.numThreads = numberOfThreads;
	}
	
	/**
	 * Adds task to the threadPoolExecutor for this instance. If the task queue limit is exceeded,
	 * the pool will be relaunched afterwards, after waiting for all queued tasks to finish.
	 * @param task task to add to the pool
	 * @throws InterruptedException if the relauch process is interrupted
	 */
	public void queueTask(Runnable task) throws InterruptedException {
		int taskCount = pool.getQueue().size();
		if(taskCount == maxTaskCount) {
			relaunchPool();
		}
		pool.execute(task);
	}
	
	/**
	 * Terminates the pool, shutting it down and waiting for it to finish all queued tasks.
	 * @throws InterruptedException if the shutdown operation is interrupted
	 */
	public void terminatePool() throws InterruptedException  {
		pool.shutdown();
    	pool.awaitTermination(maxTaskCount*secondsPerTask, TimeUnit.SECONDS);
    	if(!pool.isShutdown()) {
			throw new InterruptedException("The ThreadPoolExecutor was not shutdown after an await Termination call");
		}
	}
	
	
	/**
	 * Shuts down the pool, waiting for all queued tasks to finish. Then, it creates a new one.
	 * @throws InterruptedException if the shutdown operation is interrupted or if the pool was not shutdown correctly.
	 */
	private void relaunchPool() throws InterruptedException {
		this.terminatePool();
    	//Create new pool
		pool = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
	}
}
