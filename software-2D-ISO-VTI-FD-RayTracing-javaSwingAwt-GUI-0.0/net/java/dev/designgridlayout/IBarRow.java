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
 * Any row created by {@code DesignGridLayout.row().bar()} implements this
 * interface. Through this interface, you can add components to the current row.
 * <p/>
 * All added components share the same width (the maximum of all components
 * preferred widths) and can be aligned on the left, center or right depending on
 * which method is called to add them, or the {@link Tag} they are assigned when
 * added.
 * <p/>
 * The actual position of added components depends on the current platform;
 * standard buttons (as determined by their assigned {@link Tag}) will always follow
 * the guidelines for the current platform.
 *
 * @author Jean-Francois Poilpret
 */
public interface IBarRow extends IHideable
{
	/**
	 * Adds components to the left of this row. Components are added left to 
	 * right, in the same order as they appear in the arguments list.
	 * <p/>
	 * This method is equivalent to:
	 * <pre>
	 * for (JComponent child: children)
	 * {
	 *     add(child, Tag.LEFT);
	 * }
	 * </pre>
	 * 
	 * @param children components to add to the left of this row
	 * @return {@code this} row (to allow chaining other methods for the current 
	 * row)
	 */
	public abstract IBarRow left(JComponent... children);
	
	/**
	 * Adds components to the center of this row. Components are added left to 
	 * right, in the same order as they appear in the arguments list.
	 * <p/>
	 * This method is equivalent to:
	 * <pre>
	 * for (JComponent child: children)
	 * {
	 *     add(child, Tag.OTHER);
	 * }
	 * </pre>
	 * 
	 * @param children components to add to the center of this row
	 * @return {@code this} row (to allow chaining other methods for the current 
	 * row)
	 */
	public abstract IBarRow center(JComponent... children);
	
	/**
	 * Adds components to the right of this row. Components are added left to 
	 * right, in the same order as they appear in the arguments list.
	 * <p/>
	 * This method is equivalent to:
	 * <pre>
	 * for (JComponent child: children)
	 * {
	 *     add(child, Tag.RIGHT);
	 * }
	 * </pre>
	 * 
	 * @param children components to add to the right of this row
	 * @return {@code this} row (to allow chaining other methods for the current 
	 * row)
	 */
	public abstract IBarRow right(JComponent... children);
	
	/**
	 * Adds a tagged component to this row. The component is placed in the row at a
	 * location determined by guidelines for the current platform, based on its
	 * {@link Tag}.
	 * 
	 * @param child the component to add to this row
	 * @param tag the tag to assign to {@code child} that will determine the
	 * actual position of {@code child} in the row
	 * @return {@code this} row (to allow chaining other methods for the current 
	 * row)
	 */
	public abstract IBarRow add(JComponent child, Tag tag);
	
	/**
	 * Adds a larger gap (see {@link javax.swing.LayoutStyle.ComponentPlacement#UNRELATED})
	 * after the component that has been added just before this method is called. Has no
	 * effect if no component was added before. Multiple consecutive calls have the same
	 * effect as only one call, DesignGridLayout will never accumulate consecutive gaps.
	 * 
	 * @return {@code this} row (to allow chaining other methods for the current 
	 * row)
	 */
	public abstract IBarRow gap();
	
	/**
	 * Makes this row independent of other non-grid rows in terms of component 
	 * width. By default, DesignGridLayout ensures that all non-grid rows in a 
	 * layout use the same width for all components of all these rows. In general,
	 * this is the behavior expected by end-users; however, there may be some 
	 * situations where this behavior is not desirable.
	 * <p/>
	 * Note that you can disable DesignGridLayout feature of consistent widths across 
	 * non-grid rows with {@link DesignGridLayout#withoutConsistentWidthAcrossNonGridRows()}.
	 * 
	 * @return {@code this} row (to allow chaining other methods for the current 
	 * row)
	 */
	public abstract IBarRow withOwnRowWidth();
}
