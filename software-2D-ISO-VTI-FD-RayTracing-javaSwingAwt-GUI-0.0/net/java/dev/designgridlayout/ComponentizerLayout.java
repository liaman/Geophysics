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

import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.LayoutManager2;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.LayoutStyle.ComponentPlacement;

import net.java.dev.designgridlayout.Componentizer.Builder;
import net.java.dev.designgridlayout.Componentizer.WidthPolicy;

final class ComponentizerLayout implements LayoutManager2, Builder
{
	ComponentizerLayout(JComponent parent)
	{
		_wrapper = new ParentWrapper<JComponent>(parent);
		_orientation = new OrientationPolicy(parent);
		parent.setLayout(this);
	}
	
	@Override public Builder withSmartVerticalResize()
	{
		_heightTester = _defaultHeightTester;
		return this;
	}

	@Override public Builder withoutSmartVerticalResize()
	{
		_heightTester = new UnitHeightGrowPolicy(_defaultHeightTester);
		return this;
	}

	@Override public Builder add(WidthPolicy width, JComponent... children)
	{
		checkNoDuplicateComponents(children);
		for (JComponent component: children)
		{
			_wrapper.checkAddedComponent(component);
		}
		switch (width)
		{
			case MIN_TO_PREF:
			_numComponentsWiderThanMin += children.length;
			break;
			
			case MIN_AND_MORE:
			_numComponentsWiderThanMin += children.length;
			_numComponentsWiderThanPref += children.length;
			break;
			
			case PREF_AND_MORE:
			_numComponentsWiderThanPref += children.length;
			break;
			
			default:
			break;
		}
		for (JComponent child: children)
		{
			_children.add(new ComponentizerItem(child, width));
			_wrapper.add(child);
		}
		return this;
	}
	
	static private void checkNoDuplicateComponents(JComponent[] components)
	{
		for (JComponent comp1: components)
		{
			int count = 0;
			for (JComponent comp2: components)
			{
				if (comp1 == comp2)
				{
					count++;
				}
			}
			if (count > 1)
			{
				throw new IllegalArgumentException("Do not add the same component twice");
			}
		}
	}
	
	@Override public Builder fixedPref(JComponent... children)
	{
		return add(WidthPolicy.PREF_FIXED, children);
	}

	@Override public Builder prefAndMore(JComponent... children)
	{
		return add(WidthPolicy.PREF_AND_MORE, children);
	}
	
	@Override public Builder minToPref(JComponent... children)
	{
		return add(WidthPolicy.MIN_TO_PREF, children);
	}
	
	@Override public Builder minAndMore(JComponent... children)
	{
		return add(WidthPolicy.MIN_AND_MORE, children);
	}
	
	@Override public JComponent component()
	{
		return _wrapper.parent();
	}

	int getBaseline()
	{
		computeAll();
		return _baseline;
	}
	
	@Override public void layoutContainer(Container parent)
	{
		_wrapper.checkParent(parent);
		
		synchronized(parent.getTreeLock())
		{
			computeAll();

			int parentWidth = parent.getSize().width;
			// Never layout components smaller than the minimum size
			parentWidth = Math.max(parentWidth, _minWidth);
			int availableWidth = parentWidth - _gap;
			int prefWidth = _prefWidth - _gap;

			// Prepare layout
			LayoutHelper helper = new LayoutHelper(
				_heightTester, parentWidth, _orientation.isRightToLeft());
			helper.setRowAvailableHeight(_wrapper.parent().getHeight());
			helper.setY(0);

			if (availableWidth < prefWidth)
			{
				// - if available width < pref width, use "min-pref" width of all components
				int minWidth = _minWidth - _gap;
				float ratio = 0.0f;
				if (availableWidth > minWidth)
				{
					// Calculate extra width for each variable width component
					int extra = availableWidth - minWidth;
					ratio = (float) extra / (prefWidth - minWidth);
				}
				layoutComponentsMinToPref(helper, ratio);
			}
			else
			{
				// - if available width > pref width, increase size only of variable components
				int extra = 0;
				int fudge = 0;
				if (_numComponentsWiderThanPref > 0)
				{
					// Calculate extra width for each variable width component
					extra = (availableWidth - prefWidth) / _numComponentsWiderThanPref;
					// Fudge is assigned to the first variable width components
					// (1 pixel per component)
					fudge = (availableWidth - prefWidth) % _numComponentsWiderThanPref;
				}
				layoutComponentsPrefAndMore(helper, extra, fudge);
			}
		}
	}
	
	private void layoutComponentsMinToPref(LayoutHelper helper, float ratio)
	{
		int nth = 0;
		int x = 0;
		
		// Perform actual layout
		for (ComponentizerItem child: _children)
		{
			int width = child.minimumWidth();
			int extra = child.preferredWidth() - width;
			if (extra > 0 && needExtraWidth(child, true))
			{
				width += extra * ratio;
			}
			helper.setSizeLocation(child.component(), x, width, _height, _baseline);
			x += width + _gaps[nth];
			nth++;
		}
	}
	
	private void layoutComponentsPrefAndMore(LayoutHelper helper, int extra, int fudge)
	{
		int nth = 0;
		int x = 0;
		int remainingFudge = fudge;
		
		// Perform actual layout
		for (ComponentizerItem child: _children)
		{
			int width = child.preferredWidth();
			if (extra > 0 && needExtraWidth(child, false))
			{
				width += extra;
				if (remainingFudge > 0)
				{
					width++;
					remainingFudge--;
				}
			}
			helper.setSizeLocation(child.component(), x, width, _height, _baseline);
			x += width + _gaps[nth];
			nth++;
		}
	}
	
	static private boolean needExtraWidth(ComponentizerItem child, boolean minToPref)
	{
		switch (child.widthPolicy())
		{
			case MIN_AND_MORE:
			return true;
			
			case MIN_TO_PREF:
			return minToPref;
			
			case PREF_AND_MORE:
			return !minToPref;
			
			case PREF_FIXED:
			default:
			return false;
		}
	}

	@Override public Dimension minimumLayoutSize(Container parent)
	{
		initSizeCalculation(parent);
		return new Dimension(_minWidth, _height);
	}

	@Override public Dimension preferredLayoutSize(Container parent)
	{
		initSizeCalculation(parent);
		return new Dimension(_prefWidth, _height);
	}
	
	@Override public Dimension maximumLayoutSize(Container parent)
	{
		initSizeCalculation(parent);
		int maxHeight = (_variableHeight ? Integer.MAX_VALUE: _height);
		return new Dimension(Integer.MAX_VALUE, maxHeight);
	}

	@Override public float getLayoutAlignmentX(Container target)
	{
		return LAYOUT_ALIGNMENT;
	}

	@Override public float getLayoutAlignmentY(Container target)
	{
		return LAYOUT_ALIGNMENT;
	}

	@Override public void invalidateLayout(Container target)
	{
		_inited = false;
	}
	
	private void initSizeCalculation(Container parent)
	{
		_wrapper.checkParent(parent);
		computeAll();
	}

	@Override public void addLayoutComponent(String name, Component comp)
	{
		// This method is never called for LayoutManager2 in fact
	}

	@Override public void addLayoutComponent(Component comp, Object constraints)
	{
		_wrapper.checkAdd();
	}

	@Override public void removeLayoutComponent(Component comp)
	{
		throw new IllegalArgumentException("Do not use this method");
	}
	
	private void computeAll()
	{
		if (!_inited)
		{
			_baseline = ComponentHelper.maxValues(_children, BaselineExtractor.INSTANCE);
			_minWidth = ComponentHelper.sumValues(_children, MinWidthExtractor.INSTANCE);
			_prefWidth = ComponentHelper.sumValues(_children, PrefWidthExtractor.INSTANCE);
			_height = ComponentHelper.maxValues(_children, PrefHeightExtractor.INSTANCE);

			ComponentGapsHelper helper = ComponentGapsHelper.instance();
			_gaps = new int[_children.size()];
			_gap = 0;
			for (int nth = 0; nth < _children.size() - 1; nth++)
			{
				JComponent left = _children.get(nth).component();
				JComponent right = _children.get(nth + 1).component();
				int gap = helper.getHorizontalGap(
					left, right, ComponentPlacement.RELATED, _wrapper.parent());
				_gaps[nth] = gap;
				_gap += gap;
			}
			_minWidth += _gap;
			_prefWidth += _gap;
			if (!_children.isEmpty())
			{
				_gaps[_children.size() - 1] = 0;
			}
			_variableHeight = false;
			for (ComponentizerItem child: _children)
			{
				if (_heightTester.canGrowHeight(child.component()))
				{
					_variableHeight = true;
					break;
				}
			}

			_inited = true;
		}
	}

	static final private float LAYOUT_ALIGNMENT = 0.5f; 
	
	static private HeightGrowPolicy _defaultHeightTester = new DefaultGrowPolicy();

	private final ParentWrapper<JComponent> _wrapper;
	private final List<ComponentizerItem> _children = new ArrayList<ComponentizerItem>();
	private HeightGrowPolicy _heightTester = _defaultHeightTester;
	private final OrientationPolicy _orientation;
	private boolean _inited = false;
	private int _baseline = 0;
	private int _height = 0;
	private int _minWidth = 0;
	private int _prefWidth = 0;
	private boolean _variableHeight = false;
	private int[] _gaps = null;
	private int _gap = 0;
	private int _numComponentsWiderThanPref = 0;
	private int _numComponentsWiderThanMin = 0;
}
