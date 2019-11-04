package ngsep.gbs;

import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import ngsep.assembly.SelfAlignment;

public class ThreadPoolManager {
	private int maxTaskCount;
	private final int numThreads;
	private ThreadPoolExecutor pool;
	
	public ThreadPoolManager(int numberOfThreads, int maxTaskCount) {
		this.pool = new ThreadPoolExecutor(numberOfThreads, numberOfThreads, Long.MAX_VALUE, TimeUnit.NANOSECONDS, new LinkedBlockingQueue<Runnable>());
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
		pool.execute(task);
		int taskCount = pool.getQueue().size();
		if(taskCount > maxTaskCount) {
			relaunchPool();
		}
	}
	
	/**
	 * Terminates the pool, shutting it down and waiting for it to finish all queued tasks.
	 * @throws InterruptedException if the shutdown operation is interrupted
	 */
	public void terminatePool() throws InterruptedException  {
		pool.shutdown();
    	pool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
	}
	
	/**
	 * Terminates the pool and also ensures that the pool reference is deleted.
	 * @throws InterruptedException if the shutdown operation is interrupted
	 */
	public void killPool() throws InterruptedException {
		this.terminatePool();
		pool = null;
	}
	
	/**
	 * Shuts down the pool, waiting for all queued tasks to finish. Then, it creates a new one.
	 * @throws InterruptedException if the shutdown operation is interrupted or if the pool was not shutdown correctly.
	 */
	private void relaunchPool() throws InterruptedException {
		this.terminatePool();
    	
    	if(!pool.isShutdown()) {
			throw new InterruptedException("The ThreadPoolExecutor was not shutdown after an await Termination call");
		}
    	
    	//Create new pool
		pool = new ThreadPoolExecutor(numThreads, numThreads, Long.MAX_VALUE, TimeUnit.NANOSECONDS, new LinkedBlockingQueue<Runnable>());
	}
}
