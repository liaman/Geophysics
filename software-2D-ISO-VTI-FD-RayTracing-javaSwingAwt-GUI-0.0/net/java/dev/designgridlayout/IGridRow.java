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

import javax.swing.JComponent;

/**
 * Any row created by {@code DesignGridLayout.row().grid()} implements this
 * interface. Through this interface, you can add components to the current row;
 * position and size of components will use the canonical grid calculated for
 * the whole layout.
 * 
 * @author Jean-Francois Poilpret
 */
public interface IGridRow extends IRow, ISubGridStarter
{
	/*
	 * (non-Javadoc)
	 * @see IRow#add(javax.swing.JComponent[])
	 */
	@Override public abstract IGridRow add(JComponent... children);
	
	/**
	 * Adds one component to this row and allows it to span several columns of
	 * the canonical grid.
	 * <p/>
	 * The order of calls match the order in which the components will be added
	 * to this row. Components are added left to right, in the same order as 
	 * they appear in the arguments list.
	 * 
	 * @param child component to add to this row; it is added to the right of
	 * the component that was added by the nearest previous call to an add
	 * method.
	 * @param span the number of columns to span (must be &gt; 0)
	 * @return {@code this} row (to allow chaining other methods for the current 
	 * row)
	 */
	public abstract IGridRow add(JComponent child, int span);

	/**
	 * Adds an empty column to the current row.
	 * 
	 * @return {@code this} row (to allow chaining other methods for the current 
	 * row)
	 */
	public abstract IGridRow empty();

	/**
	 * Adds one or more empty columns to the current row.
	 * 
	 * @param span the number of columns to span (must be &gt; 0)
	 * @return {@code this} row (to allow chaining other methods for the current 
	 * row)
	 */
	public abstract IGridRow empty(int span);

	/*
	 * (non-Javadoc)
	 * @see IRow#addMulti(javax.swing.JComponent[])
	 */
	@Override public abstract IGridRow addMulti(JComponent... children);

	/*
	 * (non-Javadoc)
	 * @see net.java.dev.designgridlayout.IRow#indent()
	 */
	@Override public abstract IGridRow indent();
	
	/*
	 * (non-Javadoc)
	 * @see net.java.dev.designgridlayout.IRow#indent(int)
	 */
	@Override public abstract IGridRow indent(int n);

	/**
	 * Adds components to this row; all components are "assembled" as one
	 * global component and span a given number of columns in the row.
	 * <p/>
	 * Note that the width of each individual component will never grow bigger
	 * than its preferred width.
	 * 
	 * @param span the number of columns to span (must be &gt; 0)
	 * @param children components to assemble and add to this row
	 * @return {@code this} row (to allow chaining other methods for the current 
	 * row)
	 * 
	 * @deprecated Use {@link #add(JComponent, int)} with {@link Componentizer} instead.
	 */
	@Deprecated
	public abstract IGridRow addMulti(int span, JComponent... children);
}
