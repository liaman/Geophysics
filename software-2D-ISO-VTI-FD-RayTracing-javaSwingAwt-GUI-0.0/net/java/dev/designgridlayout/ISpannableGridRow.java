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
 * Rows created by {@code DesignGridLayout.row().grid()} and 
 * {@code DesignGridLayout.row().grid(JLabel)} implement this interface.
 * Through this interface, in addition to methods inherited from {@link IGridRow},
 * you can specify components spanning multiple rows with the new
 * {@link #spanRow()} method.
 * 
 * @author Jean-Francois Poilpret
 */
public interface ISpannableGridRow extends IGridRow
{
	/**
	 * Specifies that the next grid column will be filled with the same
	 * component in the same column location in the row above this row.
	 * <p/>
	 * Calling {@code spanRow()} in one row behaves as follows:
	 * <ol>
	 * <li>the component on the previous row at the same column position in the 
	 * same sub-grid will span this row</li>
	 * <li>the horizontal span of this component will be the same as the column 
	 * span defined for that same component on the previous row</li>
	 * <li>the sub-grid in which {@code spanRow()} appears will have its 
	 * gridspan forced to the same value as the sub-grid in the same position of 
	 * the previous row</li>
	 * <li>if the matching component in the previous row is {@link IGridRow#empty()}
	 * or {@link IGridRow#empty(int)}, then {@code spanRow()} is also {@code empty()}
	 * or {@code empty(int)}</li>
	 * <li>if this is the first row in the layout, then {@code spanRow()} is 
	 * replaced with a marker component, used for fixing the layout</li>
	 * <li>if there is no grid row above {@code this} row, then {@code spanRow()}
	 * will be replaced with a marker component, used for fixing the layout</li>
	 * <li>if there is no matching subgrid (same position) in the above row, then 
	 * {@code spanRow()} will be replaced with a marker component, used for 
	 * fixing the layout</li>
	 * <li>if there is no matching component (same position) in the matching
	 * subgrid of the above row, then {@code spanRow()} will be replaced with a 
	 * marker component, used for fixing the layout</li>
	 * <li>if the matching subgrid in the above row does not have the same number
	 * of columns of the current subgrid in {@code this} row, then 
	 * {@code spanRow()} will be replaced with a marker component, used for 
	 * fixing the layout</li>
	 * </ol>
	 * <p/>
	 * A marker component is used to show bad calls to {@code spanRow()} which
	 * could not be avoided at compile-time. This marker is a "<b>spanRow()</b>" 
	 * label with a red background and a tooltip giving more details about the 
	 * actual problem.
	 * 
	 * @return {@code this} row (to allow chaining other methods for the current 
	 * row)
	 */
	public abstract ISpannableGridRow spanRow();

	/*
	 * (non-Javadoc)
	 * @see net.java.dev.designgridlayout.IGridRow#add(javax.swing.JComponent[])
	 */
	@Override public abstract ISpannableGridRow add(JComponent... children);
	
	/*
	 * (non-Javadoc)
	 * @see net.java.dev.designgridlayout.IGridRow#add(javax.swing.JComponent, int)
	 */
	@Override public abstract ISpannableGridRow add(JComponent child, int span);
	
	/*
	 * (non-Javadoc)
	 * @see net.java.dev.designgridlayout.IGridRow#empty()
	 */
	@Override public abstract ISpannableGridRow empty();
	
	/*
	 * (non-Javadoc)
	 * @see net.java.dev.designgridlayout.IGridRow#empty(int)
	 */
	@Override public abstract ISpannableGridRow empty(int span);
	
	/*
	 * (non-Javadoc)
	 * @see net.java.dev.designgridlayout.IGridRow#addMulti(javax.swing.JComponent[])
	 */
	@Override public abstract ISpannableGridRow addMulti(JComponent... children);
	
	/*
	 * (non-Javadoc)
	 * @see net.java.dev.designgridlayout.IGridRow#addMulti(int, javax.swing.JComponent[])
	 */
	@Override public abstract ISpannableGridRow addMulti(int span, JComponent... children);

	/*
	 * (non-Javadoc)
	 * @see net.java.dev.designgridlayout.IGridRow#indent()
	 */
	@Override public abstract ISpannableGridRow indent();
	
	/*
	 * (non-Javadoc)
	 * @see net.java.dev.designgridlayout.IGridRow#indent(int)
	 */
	@Override public abstract ISpannableGridRow indent(int n);
}
