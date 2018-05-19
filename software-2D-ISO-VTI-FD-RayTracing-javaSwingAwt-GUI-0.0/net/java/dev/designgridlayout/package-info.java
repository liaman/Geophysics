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

/**
 * {@code DesignGridLayout} is a revolutionary {@link java.awt.LayoutManager} for
 * Swing applications, based on "canonical grids", well-known in publishing.
 * <p/>
 * DesignGridLayout allows to define dialogs that are always visually balanced, 
 * without ever needing any graphical designer, thanks to its fluent and 
 * easy-to-use API, which enables you to literally "visualize" the layout in the 
 * code itself.
 * <p/>
 * Using DesignGridLayout is straightforward:
 * <pre>
 * public class MyPanel extends JPanel {
 *     public MyPanel() {
 *         DesignGridLayout layout = new DesignGridLayout(this);
 *         
 *         layout.row().grid(labelA).add(fieldA);
 *         layout.row().grid(labelB).add(fieldB);
 *         //...
 *         layout.row().bar().add(okButton, Tag.OK).add(cancelButton, Tag.CANCEL);
 *     }
 *     
 *     private JLabel labelA = new JLabel("aaa");
 *     private JTextField fieldA = new JTextField();
 *     private JLabel labelB = new JLabel("bbb");
 *     private JTextField fieldB = new JTextField();
 *     //...
 *     private JButton okButton = new JButton("OK");
 *     private JButton cancelButton = new JButton("Cancel");
 * }
 * </pre>
 * As you can see in this example, each row of components in the panel can be 
 * set through one single line of source code. Each component is added from left 
 * to right (except for the last rows which is special).
 * <p/>
 * Labels (when created with 
 * {@link net.java.dev.designgridlayout.ISubGridStarter#grid(javax.swing.JLabel)})
 * have a special treatment: they are automatically right-aligned and the 
 * label column in the panel has constant width.
 * <p/>
 * All gaps between components and between components and panel borders are
 * automatically calculated by DesignGridLayout according to the current 
 * installed LookAndFeel: no need to hardcode anything!
 * <p/>
 * DesignGridLayout offers 5 types of rows:
 * <ul>
 *  <li>Grid Row: all components added to this kind of row are visually balanced
 *  	according to other grid rows, and follow best practices applied in
 *  	publishing on how to split the panel among "columns" of components</li>
 *  <li>Center Row: all components in this row are centered, they are never sized
 *  	bigger than their preferred size (except when {@code fill()} is
 *  	used (see below))</li>
 *  <li>Left Row: all components in this row are aligned on the left, they are
 *  	never sized bigger than their preferred size (except when {@code fill()} 
 *  	is used (see below))</li>
 *  <li>Right Row: all components in this row are aligned on the right, they are
 *  	never sized bigger than their preferred size (except when {@code fill()}
 *  	is used (see below))</li>
 *  <li>Command Bar Row: this row is used for command buttons such as "OK",
 *  	"Cancel", "Help"... Each button is added with a special {@code Tag},
 *  	which DesignGridLayout uses to determine the correct position of each
 *  	button based on the current platform. In this row, thus, components are
 *  	not located in the same order as they are added</li>
 *  <li>Empty Row: this row has no component but has a fixed height (automatically
 *  	calculated to "visually separate" two groups of components). This is 
 *  	useful when you want to add some space between groups of rows.</li>
 * </ul>
 * Center, Left and Right Rows have a special "fill" option that allows their 
 * extreme component(s) (rightmost component for a Left Row, leftmost component 
 * for a Right Row, and both leftmost and rightmost components for a Center Row) 
 * fill all the extra space. This can be useful for instance when you want to 
 * visually separate groups of rows:
 * <pre>
 *     layout.row().left().fill().add(new JLabel("Group"), new JSeparator());
 * </pre>
 * <p/>
 * DesignGridLayout allows you to add empty rows with a height that is 
 * automatically calculated (depending on the current installed Look &amp; Feel),
 * in order to introduce some space in your layout (e.g. to separate different 
 * groups of logically related items):
 * <pre>
 *     layout.row().grid().add(...);
 *     layout.emptyRow();
 *     layout.row().bar().add(new JButton("OK"), Tag.OK).add(new JButton("Cancel"), Tag.CANCEL);
 * </pre>
 * In grid rows (added by calling {@code DesignGridLayout.row().grid()}), you 
 * may specify that a given component spans several columns, this way, you can 
 * ensure that fields that require longer input occupy enough space to render 
 * this input, compared with shorter fields):
 * <pre>
 *     layout.row().grid().add(new JTextField(), new JTextField());
 *     layout.row().grid().add(new JTextField(), 2).add(new JTextField());
 *     layout.row().grid().add(new JTextField()).empty();
 * </pre>
 * In this snippet, the first row has two short text fields on two columns (one
 * per field); the second row has one long text field (on two columns) and one 
 * short text field (on one column).
 * <p/>
 * For any Grid Row, the number of columns is normally defined by the number of
 * added components in the row, components being counted as many times as their
 * associated span (when explicitly specified). Note however, that you can also
 * introduce empty columns in such a row:
 * <pre>
 *     layout.row().grid().empty().add(new JTextField()).empty(2);
 * </pre>
 * This code creates a row with four columns, but only the second contains a 
 * real component.
 * 
 * <h3>Forms with multiple label columns</h3>
 * With DesignGridLayout, it is also possible to define layouts with several 
 * label columns. 
 * <p/>
 * Each such label column actually defines a new canonical grid (called 
 * <b>"sub-grid"</b>). In one grid row, you can indicate you want to start a new 
 * sub-grid by calling one of the existing 
 * {@link net.java.dev.designgridlayout.ISubGridStarter#grid()} methods, then
 * all components added to this sub-grid (by one of the several 
 * {@link net.java.dev.designgridlayout.IGridRow IGridRow#add()} methods) will 
 * be part of the sub-grid:
 * <pre>
 *     layout.row().grid(label1).add(field1, field2).grid(label3).add(field3, field4);
 *     layout.row().grid(label5).add(field5)        .grid()      .add(field6);
 * </pre>
 * In this snippet, there are two rows and two sub-grids; in the first sub-grid, 
 * there are 3 components, {@code field1}, {@code field2}, {@code field5}, in
 * addition to labels {@code label1} and {@code label5}; the second sub-grid
 * contains 3 components, {@code field3}, {@code field4}, {@code field6}, and
 * one label {@code label5}; note that the second sub-grid in the second row
 * has no label.
 * <p/>
 * To handle multiple sub-grids, special processing occurs in DesignGridLayout:
 * <ul>
 * <li>first of all, find out the total number of sub-grids in the whole layout
 * (this is defined by the maximum number of calls to {@code grid()} among all
 * rows of the layout). Note that this calculus is a bit more complex because
 * it also accounts for {@code gridspan}s specified when starting a sub-grid).</li>
 * <li>for each sub-grid, the width of its label column is always fixed 
 * (calculated based on the largest label in the column); the other components
 * in the sub-grid are organized canonically. All canonical parts of all 
 * sub-grids (thus excluding the fixed-width label column) get an equal part of
 * the available width.</li>
 * </ul>
 * Each sub-grid can be assigned an explicit {@code gridspan} to specify that it
 * will span several sub-grids of the other rows. If you don't specify any
 * {@code gridspan} when starting a new sub-grid, then this sub-grid will 
 * automatically span all space on its right unless another call to a 
 * {@code grid()} method occurs in the same row:
 * <pre>
 *     layout.row().grid(label1).add(field1, field2).grid(label3).add(field3);
 *     layout.row().grid(label4).add(field4);
 * </pre>
 * In this snippet, there are 2 rows and 2 sub-grids; the first sub-grid in the
 * second row ({@code field4}) spans the whole width of the form. If you want 
 * the same layout but don't want {@code field4} to span the second sub-grid, 
 * then you can write:
 * <pre>
 *     layout.row().grid(label1)   .add(field1, field2).grid(label3).add(field3);
 *     layout.row().grid(label4, 1).add(field4);
 * </pre>
 * If you want more details about this feature, you can look at the API 
 * documentation for {@link net.java.dev.designgridlayout.ISubGridStarter}.
 * 
 * <h3>Components spanning multiple grid rows</h3>
 * DesignGridLayout also allows you to define components that span several rows.
 * This feature is available only on grid-rows (i.e. rows created with 
 * {@code layout.row().grid(...)}).
 * The way to use it is quite straightforward, you define a first row as usual:
 * <pre>
 *     layout.row().grid(label1).add(field1).add(list);
 * </pre>
 * Then on the <b>next</b> row, you specify which grid column will "share" the
 * same component as the first row (that component will thus span row 2):
 * <pre>
 *     layout.row().grid(label2).add(field2).spanRow();
 * </pre>
 * Please note the use of {@code spanRow()} (instead of {@code add(...)}) to
 * specify that we want the matching component in the previous row to span this
 * row.
 * <p/>
 * This feature is described in further details in the API documentation of
 * {@link net.java.dev.designgridlayout.ISpannableGridRow#spanRow()}.
 * <p/>
 * There are some limitations in the use of {@code spanRow()} that, unfortunately,
 * are not detectable at compile-time. When you mistakenly call {@code spanRow()}
 * in an impossible layout situation, DesignGridLayout will not throw any
 * Exception at run-time, but will replace your calls to {@code spanRow()} by
 * a special marker component: a <b>"spanRow()"</b> label with red background 
 * and a tooltip describing the problem further.
 * 
 * @author Jason Aaron Osgood
 * @author Jean-Francois Poilpret
 * @see net.java.dev.designgridlayout.DesignGridLayout
 */
package net.java.dev.designgridlayout;
