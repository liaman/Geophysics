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
 * Utility to create a new {@link JComponent} out of several {@link JComponent}s,
 * assembled horizontally, separated with natural gaps (as provided by the current
 * look and feel), and which minimum and preferred sizes are calculated based on
 * the {@link WidthPolicy} defined by each component when it's added.
 * <p/>
 * {@link Componentizer} is particularly useful when used in conjunction with
 * {@link DesignGridLayout}, when you would like to add a set of components into
 * only one cell (one column), but it can also be used with other 
 * {@link LayoutManager}s.
 * <p/>
 * The following snippet shows creation of a component made of a {@code JTextField}
 * followed by a {@code JLabel}, used to represent a quantity followed by its measure 
 * unit:
 * <pre>
 * // Create the widget for pressure (quantity + unit) 
 * JTextField quantity = new JTextField("1080");
 * JLabel unit = new JLabel("hPa"):
 * JComponent pressure = Componentizer.create()
 *     .add(WidthPolicy.PREF_AND_MORE, quantity)
 *     .add(WidthPolicy.PREF_FIXED, unit)
 *     .component();
 * // Use the pressure widget in any layout
 * ...
 * </pre>
 * 
 * @author Jean-Francois Poilpret
 */
final public class Componentizer
{
	private Componentizer()
	{
	}

	/**
	 * Create a new {@link Builder} that will be used to build a new {@link JComponent}
	 * by assembly of several {@code JComponent}s.
	 * <p/>
	 * The returned {@link Builder} is then used, with its fluent API, to add individual
	 * {@link JComponent}s.
	 * 
	 * @return a new {@link Builder}
	 */
	static public Builder create()
	{
		return new ComponentizerLayout(new MultiComponent());
	}

	/**
	 * Fluent interface to "componentize" a set of {@link JComponent}s, returned
	 * by {@link Componentizer#create()}.
	 * 
	 * @author Jean-Francois Poilpret
	 */
	static public interface Builder
	{
		/**
		 * Enable "smart vertical resize" for components that can grow in height. Smart
		 * resize means that no row can be truncated on e.g. a {@code JTable}.
		 * 
		 * @return {@code this} instance of {@code Componentizer.Builder}, allowing for 
		 * chained calls to other methods (also known as "fluent API")
		 */
		public Builder withSmartVerticalResize();
		
		/**
		 * Disable "smart vertical resize" for components that can grow in height. Smart
		 * resize means that no row can be truncated on e.g. a {@code JTable}.
		 * 
		 * @return {@code this} instance of {@code Componentizer.Builder}, allowing for 
		 * chained calls to other methods (also known as "fluent API")
		 */
		public Builder withoutSmartVerticalResize();

		/**
		 * Add {@code children} to the component that will be created by this 
		 * {@code Builder}; each added component will use the specified {@code width}
		 * as a {@link WidthPolicy} used during resize.
		 * 
		 * @param width the {@link WidthPolicy} to use for {@code children}
		 * @param children the components to assemble into one new component
		 * @return {@code this} instance of {@code Componentizer.Builder}, allowing for 
		 * chained calls to other methods (also known as "fluent API")
		 */
		public Builder add(WidthPolicy width, JComponent... children);

		/**
		 * Add {@code children} to the component that will be created by this 
		 * {@code Builder}; each added component will always keep its preferred
		 * width.
		 * <p/>
		 * This method is equivalent to {@code add(WidthPolicy.PREF_FIXED, children)}.
		 * 
		 * @param children the components to assemble into one new component
		 * @return {@code this} instance of {@code Componentizer.Builder}, allowing for 
		 * chained calls to other methods (also known as "fluent API")
		 */
		public Builder fixedPref(JComponent... children);

		/**
		 * Add {@code children} to the component that will be created by this 
		 * {@code Builder}; each added component will have a width greater or equal
		 * to its preferred width.
		 * <p/>
		 * This method is equivalent to {@code add(WidthPolicy.PREF_AND_MORE, children)}.
		 * 
		 * @param children the components to assemble into one new component
		 * @return {@code this} instance of {@code Componentizer.Builder}, allowing for 
		 * chained calls to other methods (also known as "fluent API")
		 */
		public Builder prefAndMore(JComponent... children);

		/**
		 * Add {@code children} to the component that will be created by this 
		 * {@code Builder}; each added component will have a variable width, varying
		 * from its minimum width to its preferred width.
		 * <p/>
		 * This method is equivalent to {@code add(WidthPolicy.MIN_TO_PREF, children)}.
		 * 
		 * @param children the components to assemble into one new component
		 * @return {@code this} instance of {@code Componentizer.Builder}, allowing for 
		 * chained calls to other methods (also known as "fluent API")
		 */
		public Builder minToPref(JComponent... children);

		/**
		 * Add {@code children} to the component that will be created by this 
		 * {@code Builder}; each added component will have a width greater or equal
		 * to its minimum width, but no upper bound.
		 * <p/>
		 * This method is equivalent to {@code add(WidthPolicy.PREF_AND_MORE, children)}.
		 * 
		 * @param children the components to assemble into one new component
		 * @return {@code this} instance of {@code Componentizer.Builder}, allowing for 
		 * chained calls to other methods (also known as "fluent API")
		 */
		public Builder minAndMore(JComponent... children);
		
		/**
		 * Get the {@link JComponent} resulting from the aggregation of all components
		 * added with {@link #add(WidthPolicy, JComponent...)}.
		 * 
		 * @return the aggregated {@link JComponent}
		 */
		public JComponent component();
	}
	
	/**
	 * Policy defining how a component should be resized when its embedding component
	 * (created by {@link Componentizer}) is resized horizontally.
	 * 
	 * @author Jean-Francois Poilpret
	 */
	static public enum WidthPolicy
	{
		/**
		 * A component using this policy will always have a fixed width, equal to its
		 * preferred width.
		 */
		PREF_FIXED,
		
		/**
		 * A component using this policy will have a width that can only decrease, from
		 * its preferred width, until it reaches its minimum width.
		 */
		MIN_TO_PREF,
		
		/**
		 * A component using this policy will have a width that can only increase, from
		 * its preferred size, upwards without any maximum limit.
		 */
		PREF_AND_MORE,
		
		/**
		 * A component using this policy will have a width that can either decrease, from
		 * its preferred width, until it reaches its minimum width, or increase upwards
		 * without any maximum limit.
		 */
		MIN_AND_MORE
	}
	
	static private class MultiComponent extends JComponent
	{
		private static final long serialVersionUID = 1L;

		@Override public int getBaseline(int width, int height)
		{
			return ((ComponentizerLayout) getLayout()).getBaseline();
		}
	}
}
