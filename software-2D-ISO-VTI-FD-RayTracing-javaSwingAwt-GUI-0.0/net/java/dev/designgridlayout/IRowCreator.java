//  Copyright 2005-2011 Jason Aaron Osgood, Jean-Francois Poilpret
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package net.java.dev.designgridlayout;

/**
 * Interface returned by {@link DesignGridLayout#row()} in order to specify the
 * type of row to be actually added to the current layout (grid, centered, left-
 * or right-aligned).
 * 
 * @author Jean-Francois Poilpret
 */
public interface IRowCreator extends ISubGridStarter
{
	/**
	 * Adds the current row to a given {@code group}. Grouping rows is useful
	 * when you want to show or hide several rows at a time.
	 * <p/>
	 * It is also possible to add one row to several groups.
	 * 
	 * @param group the group to which the current row must be added
	 * @return {@code this} instance of IRowCreator, allowing for chained 
	 * calls to other methods (also known as "fluent API")
	 */
	public abstract IRowCreator group(RowGroup group);

	/**
	 * Creates a center-aligned row. Center-aligned rows are NOT canonical grids 
	 * but avoid the otherwise mandatory use of several {@code LayoutManager}s 
	 * for one single dialog.
	 * <p/>
	 * The new row is located under the previously created row, which means that
	 * each line of code using this method creates a new row and all lines can 
	 * be read as defining the visual UI of the container.
	 * 
	 * @return a new center-aligned row which API is used in a chained-way 
	 * (fluent API) to add components to the row.
	 */
	public abstract INonGridRow center();

	/**
	 * Creates a left-aligned row. Left-aligned rows are NOT canonical grids but
	 * avoid the otherwise mandatory use of several {@code LayoutManager}s for 
	 * one single dialog.
	 * <p/>
	 * The new row is located under the previously created row, which means that
	 * each line of code using this method creates a new row and all lines can 
	 * be read as defining the visual UI of the container.
	 * 
	 * @return a new left-aligned row which API is used in a chained-way (fluent 
	 * API) to add components to the row.
	 */
	public abstract INonGridRow left();

	/**
	 * Creates a right-aligned row. Right-aligned rows are NOT canonical grids 
	 * but avoid the otherwise mandatory use of several {@code LayoutManager} 
	 * for one single dialog.
	 * <p/>
	 * The new row is located under the previously created row, which means that
	 * each line of code using this method creates a new row and all lines can 
	 * be read as defining the visual UI of the container.
	 * 
	 * @return a new right-aligned row which API is used in a chained-way (fluent 
	 * API) to add components to the row.
	 */
	public abstract INonGridRow right();

	/**
	 * Creates a "command bar" row. This kind of row is NOT a canonical grid and is
	 * specially dedicated to layout rows with command buttons such as "OK", "Cancel",
	 * "Apply", "Help"...
	 * <p/>
	 * Command bar rows place standard components (marked with {@link Tag}) based on
	 * platform preferences. This means that you don't have to care about where you 
	 * must locate "OK" and "Cancel" if your application must run on both Windows and
	 * Mac OS.
	 * <p/>
	 * These rows are split into 3 sub-rows: one on the left, one
	 * in the center and one on the right.
	 * <p/>
	 * The new row is located under the previously created row, which means that
	 * each line of code using this method creates a new row and all lines can 
	 * be read as defining the visual UI of the container.
	 * 
	 * @return a new command-bar row which API is used in a chained-way (fluent 
	 * API) to add components to the row.
	 */
	public abstract IBarRow bar();
}
