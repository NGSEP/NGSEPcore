package ngsep.gbs;

import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

public class ThreadPoolManager {
	private static final int TIMEOUT_SECONDS = 30;
	
	private int maxTaskCount;
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
			terminatePool();
		}
		pool.execute(task);
	}
	
	/**
	 * Terminates the pool, shutting it down and waiting for it to finish all queued tasks.
	 * @throws InterruptedException if the shutdown operation is interrupted
	 */
	public void terminatePool() throws InterruptedException  {
		pool.shutdown();
    	pool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
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
